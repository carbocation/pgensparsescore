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
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4
1\t100\tv1\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\t./.
1\t200\tv2\tC\tT\t.\tPASS\t.\tGT\t1/1\t1/1\t0/1\t1/1
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


def read_plink_sums(path: pathlib.Path) -> list[float]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    sum_columns = [name for name in rows[0] if name.endswith("_SUM")]
    if len(sum_columns) != 1:
        raise AssertionError(f"expected one PLINK score-sum column, got {sum_columns}")
    return [float(row[sum_columns[0]]) for row in rows]


def assert_close(actual: list[float], expected: list[float], label: str) -> None:
    if len(actual) != len(expected):
        raise AssertionError(f"{label}: length mismatch {len(actual)} != {len(expected)}")
    for idx, (left, right) in enumerate(zip(actual, expected, strict=True)):
        if not math.isclose(left, right, rel_tol=1e-12, abs_tol=1e-12):
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
        frequency_prefix = tmp / "tiny-frequency"
        run(
            args.plink2,
            "--pfile",
            str(pfile),
            "--freq",
            "--out",
            str(frequency_prefix),
        )
        output = tmp / "result"
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
            "--out",
            str(output),
        )

        with gzip.open(output.with_suffix(".scores.tsv.gz"), "rt", newline="") as handle:
            output_rows = list(csv.DictReader(handle, delimiter="\t"))
        if list(output_rows[0]) != ["FID", "IID", "score1", "score2"]:
            raise AssertionError(f"unexpected output columns: {list(output_rows[0])}")
        score1 = [float(row["score1"]) for row in output_rows]
        score2 = [float(row["score2"]) for row in output_rows]
        assert_close(score1, [0.0, 2.0, 7.0, 2.0], "direct score1")
        assert_close(score2, [2.0, 1.0, 0.0, 1.0], "direct score2")
        if output.with_suffix(".work.score-major.bin").exists():
            raise AssertionError("working score-major matrix was not removed")

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
                str(frequency_prefix.with_suffix(".afreq")),
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
            assert_close(score1 if score_name == "score1" else score2, expected,
                         f"PLINK oracle {score_name}")


if __name__ == "__main__":
    main()
