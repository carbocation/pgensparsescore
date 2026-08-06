#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-only

from __future__ import annotations

import argparse
import json
import pathlib
import subprocess
import tempfile


HEADER = "SNP\tEFFECT_ALLELE\tOTHER_ALLELE\tEFFECT_ALLELE_WEIGHT\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scorer", required=True)
    args = parser.parse_args()

    with tempfile.TemporaryDirectory(prefix="pgensparsescore-progress-") as raw:
        root = pathlib.Path(raw)
        (root / "score1.tsv").write_text(
            HEADER + "v1\tG\tA\t1\n" + "v2\tT\tC\t0\n"
        )
        (root / "score2.tsv").write_text(HEADER + "v1\tA\tG\t2\n")
        (root / "manifest.tsv").write_text(
            "SCORE\tPATH\nscore1\tscore1.tsv\nscore2\tscore2.tsv\n"
        )
        (root / "variant-map.tsv").write_text(
            "SOURCE_ID\tTARGET_ID\nv1\ttarget-v1\n"
        )
        progress = root / "progress.jsonl"
        catalog = root / "scores.catalog.bin"
        completed = subprocess.run(
            [
                args.scorer,
                "compile",
                "--manifest",
                str(root / "manifest.tsv"),
                "--variant-map",
                str(root / "variant-map.tsv"),
                "--progress-jsonl",
                str(progress),
                "--progress-interval-seconds",
                "1",
                "--out",
                str(catalog),
            ],
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        if completed.returncode:
            raise RuntimeError(completed.stdout)

        events = [json.loads(line) for line in progress.read_text().splitlines()]
        if [event["sequence"] for event in events] != list(range(len(events))):
            raise AssertionError("progress event sequence is not contiguous")
        phases = {event["phase"] for event in events}
        required = {
            "start",
            "variant_filter_loaded",
            "manifest_loaded",
            "sort_weights",
            "sort_variants",
            "weights_ready",
            "serialize",
            "serialization_complete",
            "complete",
        }
        if not required.issubset(phases):
            raise AssertionError(f"missing progress phases: {required - phases}")
        final = events[-1]
        if final["phase"] != "complete":
            raise AssertionError("final progress event is not complete")
        for key, expected in {"scores": 2, "variants": 1, "weights": 2}.items():
            if final[key] != expected:
                raise AssertionError(f"unexpected final {key}: {final[key]}")
        for event in events:
            for key in (
                "schema_version",
                "timestamp_unix_ms",
                "elapsed_ms",
                "operation",
                "phase",
                "rss_bytes",
                "peak_rss_bytes",
            ):
                if key not in event:
                    raise AssertionError(f"progress event lacks {key}")

        invalid = subprocess.run(
            [
                args.scorer,
                "compile",
                "--manifest",
                str(root / "manifest.tsv"),
                "--progress-interval-seconds",
                "0",
                "--out",
                str(root / "invalid.catalog.bin"),
            ],
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        if invalid.returncode == 0 or "requires --progress-jsonl" not in invalid.stdout:
            raise AssertionError(f"invalid progress options were accepted:\n{invalid.stdout}")


if __name__ == "__main__":
    main()
