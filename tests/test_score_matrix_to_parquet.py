#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-only

from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pyarrow.parquet as pq


SCRIPT = Path(__file__).parents[1] / "tools" / "score_matrix_to_parquet.py"
SPEC = importlib.util.spec_from_file_location("score_matrix_to_parquet", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class ScoreMatrixToParquetTest(unittest.TestCase):
    def test_conversion_preserves_sample_rows_and_named_scores(self) -> None:
        with tempfile.TemporaryDirectory() as raw_directory:
            directory = Path(raw_directory)
            matrix = np.array(
                [[1.25, 2.5, 3.75], [-4.0, 0.5, 8.0]], dtype="<f8"
            )
            matrix.tofile(directory / "scores.f64le")
            (directory / "samples.tsv").write_text(
                "FID\tIID\nfamily1\tsample1\nfamily2\tsample2\nfamily3\tsample3\n",
                encoding="utf-8",
            )
            with (directory / "score-metadata.tsv").open(
                "w", encoding="utf-8", newline=""
            ) as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["INDEX", "SCORE"])
                writer.writerow([0, "pgsc__PGS000001"])
                writer.writerow([1, "pgsc__PGS000002"])
            metadata = {
                "format": "pgensparsescore-score-major-f64-v1",
                "path": "scores.f64le",
                "samples_path": "samples.tsv",
                "dtype": "float64",
                "byte_order": "little-endian",
                "matrix_layout": "score-major",
                "sample_id_columns": ["FID", "IID"],
                "sample_rows": 3,
                "score_columns": 2,
            }
            (directory / "scores.json").write_text(
                json.dumps(metadata), encoding="utf-8"
            )

            result = MODULE.convert(
                directory / "scores.json",
                directory / "score-metadata.tsv",
                directory / "scores.parquet",
                batch_rows=2,
                compression="zstd",
                compression_level=3,
            )
            table = pq.read_table(directory / "scores.parquet")
            self.assertEqual(result["sample_rows"], 3)
            self.assertEqual(
                table.column_names,
                ["FID", "IID", "pgsc__PGS000001", "pgsc__PGS000002"],
            )
            self.assertEqual(table.column("IID").to_pylist(), ["sample1", "sample2", "sample3"])
            self.assertEqual(table.column("pgsc__PGS000001").to_pylist(), [1.25, 2.5, 3.75])
            self.assertEqual(table.column("pgsc__PGS000002").to_pylist(), [-4.0, 0.5, 8.0])
            self.assertEqual(pq.ParquetFile(directory / "scores.parquet").num_row_groups, 2)


if __name__ == "__main__":
    unittest.main()
