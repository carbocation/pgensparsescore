#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-only

"""Exercise the CLI on a real PGEN and compare every score to PLINK 2."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import pathlib
import subprocess
import tempfile


VCF = """\
##fileformat=VCFv4.2
##contig=<ID=1,length=1000000>
##contig=<ID=2,length=1000000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4
1\t100\tv1\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\t./.
2\t200\tv2\tC\tT\t.\tPASS\t.\tGT\t1/1\t1/1\t0/1\t1/1
"""

WEIGHT_HEADER = "SNP\tEFFECT_ALLELE\tOTHER_ALLELE\tEFFECT_ALLELE_WEIGHT\n"


def run(*args: str) -> None:
    completed = subprocess.run(
        args, check=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True
    )
    if completed.returncode:
        raise RuntimeError(
            f"command failed with exit {completed.returncode}: {' '.join(args)}\n"
            f"{completed.stdout}"
        )


def run_expect_failure(*args: str, message: str) -> None:
    completed = subprocess.run(
        args, check=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True
    )
    if completed.returncode == 0:
        raise AssertionError(f"command unexpectedly succeeded: {' '.join(args)}")
    if message not in completed.stdout:
        raise AssertionError(
            f"failed command did not report {message!r}: {' '.join(args)}\n"
            f"{completed.stdout}"
        )


def read_plink_sums(path: pathlib.Path) -> list[float]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    sum_columns = [name for name in rows[0] if name.endswith("_SUM")]
    if len(sum_columns) != 1:
        raise AssertionError(f"expected one PLINK score-sum column, got {sum_columns}")
    return [float(row[sum_columns[0]]) for row in rows]


def assert_close(
    actual: list[float], expected: list[float], label: str, *, abs_tol: float = 1e-12
) -> None:
    if len(actual) != len(expected):
        raise AssertionError(f"{label}: length mismatch {len(actual)} != {len(expected)}")
    for idx, (left, right) in enumerate(zip(actual, expected, strict=True)):
        if not math.isclose(left, right, rel_tol=1e-12, abs_tol=abs_tol):
            raise AssertionError(f"{label}[{idx}]: {left} != {right}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--plink2", required=True)
    parser.add_argument("--scorer", required=True)
    args = parser.parse_args()

    with tempfile.TemporaryDirectory(prefix="pgensparsescore-integration-") as raw_tmp:
        tmp = pathlib.Path(raw_tmp)
        (tmp / "tiny.vcf").write_text(VCF)
        (tmp / "score1.tsv").write_text(
            WEIGHT_HEADER + "v1\tG\tA\t2\n" + "v2\tC\tT\t3\n"
        )
        (tmp / "score2.tsv").write_text(WEIGHT_HEADER + "v1\tA\tG\t1\n")
        (tmp / "manifest.tsv").write_text(
            "SCORE\tPATH\nscore1\tscore1.tsv\nscore2\tscore2.tsv\n"
        )

        pfile = tmp / "tiny"
        run(
            args.plink2,
            "--vcf",
            str(tmp / "tiny.vcf"),
            "--make-pgen",
            "--out",
            str(pfile),
        )
        frequency_path = tmp / "reference.acount"
        frequency_path.write_text(
            "#CHROM\tID\tREF\tALT\tREF_CT\tALT_CTS\tOBS_CT\n"
            "1\tv1\tA\tG\t5\t1\t6\n"
            "2\tv2\tC\tT\t0\t8\t8\n"
        )
        single_output = tmp / "single-result"
        run(
            args.scorer,
            "--pgen",
            str(pfile.with_suffix(".pgen")),
            "--pvar",
            str(pfile.with_suffix(".pvar")),
            "--psam",
            str(pfile.with_suffix(".psam")),
            "--manifest",
            str(tmp / "manifest.tsv"),
            "--read-freq",
            str(frequency_path),
            "--error-on-missing-freq",
            "--out",
            str(single_output),
        )

        shard_prefixes = []
        for chromosome in ("1", "2"):
            shard = tmp / f"tiny-chr{chromosome}"
            run(
                args.plink2,
                "--pfile",
                str(pfile),
                "--chr",
                chromosome,
                "--make-pgen",
                "--out",
                str(shard),
            )
            shard_prefixes.append(shard)
        (tmp / "pfiles.tsv").write_text(
            "PGEN\tPVAR\tPSAM\n"
            + "".join(
                f"{shard.name}.pgen\t{shard.name}.pvar\t{shard.name}.psam\n"
                for shard in shard_prefixes
            )
        )
        multi_output = tmp / "multi-result"
        run(
            args.scorer,
            "--pfile-list",
            str(tmp / "pfiles.tsv"),
            "--manifest",
            str(tmp / "manifest.tsv"),
            "--read-freq",
            str(frequency_path),
            "--error-on-missing-freq",
            "--out",
            str(multi_output),
        )

        score_results = {}
        for mode, output in (("single", single_output), ("multi", multi_output)):
            with gzip.open(
                output.with_suffix(".scores.tsv.gz"), "rt", newline=""
            ) as handle:
                output_rows = list(csv.DictReader(handle, delimiter="\t"))
            if list(output_rows[0]) != ["IID", "score1", "score2"]:
                raise AssertionError(
                    f"unexpected {mode} output columns: {list(output_rows[0])}"
                )
            score1 = [float(row["score1"]) for row in output_rows]
            score2 = [float(row["score2"]) for row in output_rows]
            assert_close(score1, [0.0, 2.0, 7.0, 2.0 / 3.0], f"{mode} score1")
            assert_close(score2, [2.0, 1.0, 0.0, 5.0 / 3.0], f"{mode} score2")
            score_results[mode] = (score1, score2)
        assert_close(score_results["multi"][0], score_results["single"][0],
                     "multi-PGEN score1")
        assert_close(score_results["multi"][1], score_results["single"][1],
                     "multi-PGEN score2")
        if multi_output.with_suffix(".work.score-major.bin").exists():
            raise AssertionError("working score-major matrix was not removed")

        missing_frequency = tmp / "missing.acount"
        missing_frequency.write_text(
            "#CHROM\tID\tREF\tALT\tREF_CT\tALT_CTS\tOBS_CT\n"
            "1\tv1\tA\tG\t5\t1\t6\n"
        )
        run_expect_failure(
            args.scorer,
            "--pgen",
            str(pfile.with_suffix(".pgen")),
            "--pvar",
            str(pfile.with_suffix(".pvar")),
            "--psam",
            str(pfile.with_suffix(".psam")),
            "--manifest",
            str(tmp / "manifest.tsv"),
            "--read-freq",
            str(missing_frequency),
            "--error-on-missing-freq",
            "--out",
            str(tmp / "missing-result"),
            message="frequency file has no row for scored variant v2",
        )
        if (tmp / "missing-result.work.score-major.bin").exists():
            raise AssertionError("failed run left its working score matrix behind")

        mismatched_frequency = tmp / "mismatched.acount"
        mismatched_frequency.write_text(
            "#CHROM\tID\tREF\tALT\tREF_CT\tALT_CTS\tOBS_CT\n"
            "1\tv1\tG\tA\t5\t1\t6\n"
            "2\tv2\tC\tT\t0\t8\t8\n"
        )
        run_expect_failure(
            args.scorer,
            "--pgen",
            str(pfile.with_suffix(".pgen")),
            "--pvar",
            str(pfile.with_suffix(".pvar")),
            "--psam",
            str(pfile.with_suffix(".psam")),
            "--manifest",
            str(tmp / "manifest.tsv"),
            "--read-freq",
            str(mismatched_frequency),
            "--error-on-missing-freq",
            "--out",
            str(tmp / "mismatched-result"),
            message="frequency alleles disagree with PVAR for v1",
        )

        for score_name in ("score1", "score2"):
            plink_weights = tmp / f"{score_name}.plink.tsv"
            source_rows = (tmp / f"{score_name}.tsv").read_text().splitlines()[1:]
            plink_weights.write_text(
                "ID\tA1\tWEIGHT\n"
                + "\n".join(
                    "\t".join((row.split("\t")[0], row.split("\t")[1], row.split("\t")[3]))
                    for row in source_rows
                )
                + "\n"
            )
            oracle = tmp / f"oracle-{score_name}"
            run(
                args.plink2,
                "--pfile",
                str(pfile),
                "--read-freq",
                str(frequency_path),
                "--score",
                str(plink_weights),
                "1",
                "2",
                "3",
                "header-read",
                "cols=scoresums",
                "--out",
                str(oracle),
            )
            expected = read_plink_sums(oracle.with_suffix(".sscore"))
            assert_close(score_results["multi"][0 if score_name == "score1" else 1],
                         expected,
                         f"PLINK oracle {score_name}",
                         abs_tol=5e-6)


if __name__ == "__main__":
    main()
