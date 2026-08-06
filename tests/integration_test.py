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
import struct
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
    parser.add_argument("--dense-kernel", choices=("auto", "direct", "onemkl"),
                        default="auto")
    args = parser.parse_args()
    dense_kernel_args = ["--dense-kernel", args.dense_kernel]

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

        projected_frequency = tmp / "projected-reference.acount"
        projected_frequency.write_text(
            "#CHROM\tID\tREF\tALT\tREF_CT\tALT_CTS\tOBS_CT\n"
            "1\tv1\tA\tG\t5\t1\t6\n"
        )
        (tmp / "projected-variants.tsv").write_text(
            "SOURCE_ID\tTARGET_ID\tREF\tALT\n"
            "v1\tv1\tA\tG\n"
            "v2\tv2\tC\tT\n"
            "missing\tmissing\tG\tA\n"
        )
        projected_index = tmp / "projected.index.bin"
        run(
            args.scorer,
            "build-variant-index",
            "--variant-list",
            str(tmp / "projected-variants.tsv"),
            "--out",
            str(projected_index),
        )
        projected_support = tmp / "projected.support.bin"
        run(
            args.scorer,
            "build-support-index",
            "--variant-index",
            str(projected_index),
            "--pvar",
            str(pfile.with_suffix(".pvar")),
            "--read-freq",
            str(projected_frequency),
            "--out",
            str(projected_support),
        )
        support_metadata = json.loads(
            projected_support.with_suffix(".bin.json").read_text()
        )
        if (
            support_metadata["usable_variants"] != 1
            or support_metadata["missing_frequencies"] != 1
            or support_metadata["missing_variants"] != 1
        ):
            raise AssertionError("projected support counts are wrong")

        (tmp / "projected-score1.tsv").write_text(
            WEIGHT_HEADER
            + "v1\tG\tA\t2\n"
            + "v2\tC\tT\t3\n"
            + "missing\tA\tG\t5\n"
        )
        (tmp / "projected-score2.tsv").write_text(
            WEIGHT_HEADER + "v1\tA\tG\t1\n"
        )
        (tmp / "projected-manifest.tsv").write_text(
            "SCORE_ID\tCOLUMN_NAME\tPATH\n"
            "source:score1\tscore1\tprojected-score1.tsv\n"
            "source:score2\tscore2\tprojected-score2.tsv\n"
        )
        projected_fragment = tmp / "projected.fragment.bin"
        run(
            args.scorer,
            "compile-fragment",
            "--manifest",
            str(tmp / "projected-manifest.tsv"),
            "--variant-index",
            str(projected_index),
            "--support-index",
            str(projected_support),
            "--out",
            str(projected_fragment),
        )
        fragment_metadata = json.loads(
            projected_fragment.with_suffix(".bin.json").read_text()
        )
        if (
            fragment_metadata["catalog_weights"] != 4
            or fragment_metadata["weights"] != 2
            or fragment_metadata["missing_variant_weights"] != 1
            or fragment_metadata["missing_frequency_weights"] != 1
        ):
            raise AssertionError("projected fragment counts are wrong")
        (tmp / "projected-fragments.tsv").write_text(
            "FRAGMENT\nprojected.fragment.bin\n"
        )
        projected_output = tmp / "projected-result"
        run(
            args.scorer,
            "--pgen",
            str(pfile.with_suffix(".pgen")),
            "--pvar",
            str(pfile.with_suffix(".pvar")),
            "--psam",
            str(pfile.with_suffix(".psam")),
            "--variant-index",
            str(projected_index),
            "--support-index",
            str(projected_support),
            "--fragment-list",
            str(tmp / "projected-fragments.tsv"),
            "--read-freq",
            str(projected_frequency),
            "--missing-freq",
            "omit",
            *dense_kernel_args,
            "--out",
            str(projected_output),
        )
        with gzip.open(
            projected_output.with_suffix(".scores.tsv.gz"), "rt", newline=""
        ) as handle:
            projected_rows = list(csv.DictReader(handle, delimiter="\t"))
        assert_close(
            [float(row["score1"]) for row in projected_rows],
            [0.0, 2.0, 4.0, 2.0 / 3.0],
            "projected score1",
        )
        assert_close(
            [float(row["score2"]) for row in projected_rows],
            score_results["single"][1],
            "projected score2",
        )
        projected_binary_output = tmp / "projected-binary-result"
        run(
            args.scorer,
            "--pgen",
            str(pfile.with_suffix(".pgen")),
            "--pvar",
            str(pfile.with_suffix(".pvar")),
            "--psam",
            str(pfile.with_suffix(".psam")),
            "--variant-index",
            str(projected_index),
            "--support-index",
            str(projected_support),
            "--fragment-list",
            str(tmp / "projected-fragments.tsv"),
            "--read-freq",
            str(projected_frequency),
            "--missing-freq",
            "omit",
            *dense_kernel_args,
            "--output-format",
            "score-major-bin",
            "--out",
            str(projected_binary_output),
        )
        binary_values = struct.unpack(
            "<8d",
            projected_binary_output.with_suffix(".scores.f64le").read_bytes(),
        )
        assert_close(
            list(binary_values[:4]),
            [0.0, 2.0, 4.0, 2.0 / 3.0],
            "binary projected score1",
        )
        assert_close(
            list(binary_values[4:]),
            score_results["single"][1],
            "binary projected score2",
        )
        if projected_binary_output.with_suffix(".samples.tsv").read_text() != (
            "IID\ns1\ns2\ns3\ns4\n"
        ):
            raise AssertionError("binary sample table is wrong")
        binary_metadata = json.loads(
            projected_binary_output.with_suffix(".json").read_text()
        )
        if (
            binary_metadata["format"]
            != "pgensparsescore-score-major-f64-v1"
            or binary_metadata["path"] != "projected-binary-result.scores.f64le"
            or binary_metadata["samples_path"]
            != "projected-binary-result.samples.tsv"
            or binary_metadata["matrix_layout"] != "score-major"
        ):
            raise AssertionError("binary score metadata is wrong")
        if projected_binary_output.with_suffix(".scores.tsv.gz").exists():
            raise AssertionError("binary scoring unexpectedly wrote a wide TSV")
        if projected_binary_output.with_suffix(".scores.f64le.tmp").exists():
            raise AssertionError("binary scoring left its temporary matrix behind")
        with projected_output.with_suffix(".score-metadata.tsv").open(
            newline=""
        ) as handle:
            projected_qc = {
                row["SCORE"]: row for row in csv.DictReader(handle, delimiter="\t")
            }
        score1_qc = projected_qc["score1"]
        projected_metadata = json.loads(
            projected_output.with_suffix(".json").read_text()
        )
        if projected_metadata["dense_scoring_kernel_requested"] != args.dense_kernel:
            raise AssertionError("selected dense scoring kernel is absent from metadata")
        expected_kernel_used = "onemkl" if args.dense_kernel == "onemkl" else "direct"
        if projected_metadata["dense_scoring_kernel_used"] != expected_kernel_used:
            raise AssertionError("actual dense scoring kernel is absent from metadata")
        if args.dense_kernel == "onemkl":
            if projected_metadata["onemkl_threads"] < 1:
                raise AssertionError("oneMKL thread count is absent from metadata")
        elif projected_metadata["onemkl_threads"] != 0:
            raise AssertionError("direct scoring reports active oneMKL threads")
        if (
            score1_qc["CATALOG_WEIGHTS"] != "3"
            or score1_qc["MATCHED_WEIGHTS"] != "2"
            or score1_qc["MISSING_VARIANTS"] != "1"
            or score1_qc["MISSING_FREQUENCIES"] != "1"
            or score1_qc["SCORED_WEIGHTS"] != "1"
            or score1_qc["ALT_EFFECTS"] != "1"
            or score1_qc["REF_EFFECTS"] != "1"
        ):
            raise AssertionError(f"projected score QC is wrong: {score1_qc}")

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
        if (
            compiled_metadata["sparse_weight_edges"]
            + compiled_metadata["dense_weight_edges"]
            != compiled_metadata["weight_edges"]
        ):
            raise AssertionError("sparse/dense weight-edge accounting is wrong")
        if (
            compiled_metadata["sparse_score_updates"]
            + compiled_metadata["dense_score_updates"]
            == 0
        ):
            raise AssertionError("score-update accounting is absent")
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
