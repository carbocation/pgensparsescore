#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-only

"""Exercise the CLI on a real PGEN and compare every score to PLINK 2."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
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

        (tmp / "mapped-score1.tsv").write_text(
            WEIGHT_HEADER + "source-v1\tG\tA\t2\n" + "source-v2\tC\tT\t3\n"
        )
        (tmp / "mapped-score2.tsv").write_text(
            WEIGHT_HEADER + "source-v1\tA\tG\t1\n"
        )
        (tmp / "mapped-manifest.tsv").write_text(
            "SCORE\tPATH\nscore1\tmapped-score1.tsv\nscore2\tmapped-score2.tsv\n"
        )
        (tmp / "variant-map.tsv").write_text(
            "SOURCE_ID\tTARGET_ID\nsource-v1\tv1\nsource-v2\tv2\n"
        )
        mapped_output = tmp / "mapped-result"
        run(
            args.scorer,
            "--pfile-list",
            str(tmp / "pfiles.tsv"),
            "--manifest",
            str(tmp / "mapped-manifest.tsv"),
            "--variant-map",
            str(tmp / "variant-map.tsv"),
            "--read-freq",
            str(frequency_path),
            "--error-on-missing-freq",
            "--out",
            str(mapped_output),
        )
        with gzip.open(
            mapped_output.with_suffix(".scores.tsv.gz"), "rt", newline=""
        ) as handle:
            mapped_rows = list(csv.DictReader(handle, delimiter="\t"))
        assert_close(
            [float(row["score1"]) for row in mapped_rows],
            score_results["multi"][0],
            "variant-map score1",
        )
        assert_close(
            [float(row["score2"]) for row in mapped_rows],
            score_results["multi"][1],
            "variant-map score2",
        )
        mapped_metadata = json.loads(mapped_output.with_suffix(".json").read_text())
        if mapped_metadata["variant_mapping_rows"] != 2:
            raise AssertionError("variant mapping row count is absent from metadata")

        compiled_catalog = tmp / "mapped.catalog.bin"
        run(
            args.scorer,
            "compile",
            "--manifest",
            str(tmp / "mapped-manifest.tsv"),
            "--variant-map",
            str(tmp / "variant-map.tsv"),
            "--out",
            str(compiled_catalog),
        )
        compiled_output = tmp / "compiled-result"
        compiled_progress = tmp / "compiled-progress.jsonl"
        run(
            args.scorer,
            "--pfile-list",
            str(tmp / "pfiles.tsv"),
            "--compiled-catalog",
            str(compiled_catalog),
            "--variant-map",
            str(tmp / "variant-map.tsv"),
            "--read-freq",
            str(frequency_path),
            "--missing-freq",
            "error",
            "--progress-jsonl",
            str(compiled_progress),
            "--progress-interval-seconds",
            "1",
            "--out",
            str(compiled_output),
        )
        with gzip.open(
            compiled_output.with_suffix(".scores.tsv.gz"), "rt", newline=""
        ) as handle:
            compiled_rows = list(csv.DictReader(handle, delimiter="\t"))
        assert_close(
            [float(row["score1"]) for row in compiled_rows],
            score_results["multi"][0],
            "compiled-catalog score1",
        )
        assert_close(
            [float(row["score2"]) for row in compiled_rows],
            score_results["multi"][1],
            "compiled-catalog score2",
        )
        compiled_metadata = json.loads(
            compiled_output.with_suffix(".json").read_text()
        )
        if compiled_metadata["catalog_source"] != "compiled":
            raise AssertionError("compiled catalog source is absent from metadata")
        if compiled_metadata["pvar_variants_loaded"] != 2:
            raise AssertionError("compiled catalog did not filter PVAR metadata")
        if compiled_metadata["pgen_storage_modes"] != ["standard", "standard"]:
            raise AssertionError("ordinary PGEN storage modes are absent from metadata")
        progress_events = [
            json.loads(line) for line in compiled_progress.read_text().splitlines()
        ]
        progress_phases = {event["phase"] for event in progress_events}
        required_progress_phases = {
            "start",
            "catalog_loaded",
            "variant_map_loaded",
            "samples_loaded",
            "pvar_loaded",
            "catalog_materialized",
            "frequencies_loaded",
            "working_matrix_ready",
            "pgen_start",
            "pgen_complete",
            "pgen_scored",
            "write_scores",
            "scores_written",
            "complete",
        }
        if not required_progress_phases.issubset(progress_phases):
            raise AssertionError(
                "scoring progress lacks phases: "
                f"{required_progress_phases - progress_phases}"
            )
        final_progress = progress_events[-1]
        if (
            final_progress["phase"] != "complete"
            or final_progress["sample_rows"] != 4
            or final_progress["score_columns"] != 2
            or final_progress["pgen_inputs"] != 2
        ):
            raise AssertionError(f"invalid final scoring progress: {final_progress}")

        direct_catalog = tmp / "direct.catalog.bin"
        (tmp / "direct-score.tsv").write_text(
            WEIGHT_HEADER + "v2\tC\tT\t3\n"
        )
        (tmp / "direct-manifest.tsv").write_text(
            "SCORE\tPATH\ndirect\tdirect-score.tsv\n"
        )
        run(
            args.scorer,
            "compile",
            "--manifest",
            str(tmp / "direct-manifest.tsv"),
            "--out",
            str(direct_catalog),
        )
        direct_output = tmp / "direct-result"
        run(
            args.scorer,
            "--pgen",
            str(pfile.with_suffix(".pgen")),
            "--pvar",
            str(pfile.with_suffix(".pvar")),
            "--psam",
            str(pfile.with_suffix(".psam")),
            "--compiled-catalog",
            str(direct_catalog),
            "--read-freq",
            str(frequency_path),
            "--error-on-missing-freq",
            "--out",
            str(direct_output),
        )
        direct_metadata = json.loads(
            direct_output.with_suffix(".json").read_text()
        )
        if direct_metadata["pvar_variants_loaded"] != 1:
            raise AssertionError(
                "compiled catalog without a variant map did not filter PVAR metadata"
            )

        (tmp / "late-score.tsv").write_text(
            WEIGHT_HEADER + "source-v2\tC\tT\t3\n"
        )
        (tmp / "late-manifest.tsv").write_text(
            "SCORE\tPATH\nlate\tlate-score.tsv\n"
        )
        (tmp / "late-map.tsv").write_text(
            "SOURCE_ID\tTARGET_ID\nsource-v2\tv2\n"
        )
        late_output = tmp / "late-result"
        run(
            args.scorer,
            "--pgen",
            str(pfile.with_suffix(".pgen")),
            "--pvar",
            str(pfile.with_suffix(".pvar")),
            "--psam",
            str(pfile.with_suffix(".psam")),
            "--manifest",
            str(tmp / "late-manifest.tsv"),
            "--variant-map",
            str(tmp / "late-map.tsv"),
            "--read-freq",
            str(frequency_path),
            "--error-on-missing-freq",
            "--out",
            str(late_output),
        )
        with gzip.open(
            late_output.with_suffix(".scores.tsv.gz"), "rt", newline=""
        ) as handle:
            late_rows = list(csv.DictReader(handle, delimiter="\t"))
        assert_close(
            [float(row["late"]) for row in late_rows],
            [0.0, 0.0, 3.0, 0.0],
            "filtered PVAR original PGEN index",
        )

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

        omitted_output = tmp / "omitted-result"
        run(
            args.scorer,
            "--pgen",
            str(pfile.with_suffix(".pgen")),
            "--pvar",
            str(pfile.with_suffix(".pvar")),
            "--psam",
            str(pfile.with_suffix(".psam")),
            "--compiled-catalog",
            str(compiled_catalog),
            "--variant-map",
            str(tmp / "variant-map.tsv"),
            "--read-freq",
            str(missing_frequency),
            "--missing-freq",
            "omit",
            "--out",
            str(omitted_output),
        )
        with gzip.open(
            omitted_output.with_suffix(".scores.tsv.gz"), "rt", newline=""
        ) as handle:
            omitted_rows = list(csv.DictReader(handle, delimiter="\t"))
        assert_close(
            [float(row["score1"]) for row in omitted_rows],
            [0.0, 2.0, 4.0, 2.0 / 3.0],
            "missing-frequency omission score1",
        )
        assert_close(
            [float(row["score2"]) for row in omitted_rows],
            score_results["single"][1],
            "missing-frequency omission score2",
        )
        omitted_metadata = json.loads(omitted_output.with_suffix(".json").read_text())
        if omitted_metadata["omitted_frequency_variants"] != 1:
            raise AssertionError("omitted frequency variant count is wrong")
        with omitted_output.with_suffix(".score-metadata.tsv").open(newline="") as handle:
            omitted_score_metadata = {
                row["SCORE"]: row for row in csv.DictReader(handle, delimiter="\t")
            }
        if omitted_score_metadata["score1"]["MISSING_FREQUENCIES"] != "1":
            raise AssertionError("per-score missing frequency count is wrong")

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
