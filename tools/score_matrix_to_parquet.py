#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-only

"""Convert pgensparsescore's score-major float64 output to wide Parquet."""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
from typing import Any

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq


EXPECTED_FORMAT = "pgensparsescore-score-major-f64-v1"


def read_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"{path}: expected a JSON object")
    return value


def resolve_sibling(metadata_path: Path, value: object, field: str) -> Path:
    if not isinstance(value, str) or not value:
        raise ValueError(f"{metadata_path}: missing {field}")
    path = Path(value)
    return path if path.is_absolute() else metadata_path.parent / path


def read_samples(path: Path, expected_columns: list[str]) -> dict[str, list[str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != expected_columns:
            raise ValueError(
                f"{path}: sample columns {reader.fieldnames!r} do not match "
                f"{expected_columns!r}"
            )
        result = {name: [] for name in expected_columns}
        for row in reader:
            for name in expected_columns:
                result[name].append(row[name])
    return result


def read_score_columns(path: Path, expected_count: int) -> list[str]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"INDEX", "SCORE"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path}: missing {', '.join(sorted(missing))}")
        rows = list(reader)
    if len(rows) != expected_count:
        raise ValueError(
            f"{path}: found {len(rows)} scores; expected {expected_count}"
        )
    columns: list[str] = []
    for index, row in enumerate(rows):
        if int(row["INDEX"]) != index or not row["SCORE"]:
            raise ValueError(f"{path}: invalid score row {index + 2}")
        columns.append(row["SCORE"])
    if len(set(columns)) != len(columns):
        raise ValueError(f"{path}: duplicate score column")
    return columns


def convert(
    metadata_path: Path,
    score_metadata_path: Path,
    output_path: Path,
    batch_rows: int,
    compression: str,
    compression_level: int | None,
) -> dict[str, Any]:
    metadata = read_json(metadata_path)
    if metadata.get("format") != EXPECTED_FORMAT:
        raise ValueError(f"{metadata_path}: expected format {EXPECTED_FORMAT}")
    if metadata.get("dtype") != "float64" or metadata.get("byte_order") != "little-endian":
        raise ValueError(f"{metadata_path}: unsupported matrix dtype or byte order")
    if metadata.get("matrix_layout") != "score-major":
        raise ValueError(f"{metadata_path}: unsupported matrix layout")

    sample_count = int(metadata["sample_rows"])
    score_count = int(metadata["score_columns"])
    if sample_count < 1 or score_count < 1 or batch_rows < 1:
        raise ValueError("matrix dimensions and --batch-rows must be positive")
    sample_id_columns = metadata.get("sample_id_columns")
    if not isinstance(sample_id_columns, list) or not all(
        isinstance(value, str) and value for value in sample_id_columns
    ):
        raise ValueError(f"{metadata_path}: invalid sample_id_columns")

    matrix_path = resolve_sibling(metadata_path, metadata.get("path"), "path")
    sample_path = resolve_sibling(
        metadata_path, metadata.get("samples_path"), "samples_path"
    )
    expected_bytes = sample_count * score_count * np.dtype("<f8").itemsize
    actual_bytes = matrix_path.stat().st_size
    if actual_bytes != expected_bytes:
        raise ValueError(
            f"{matrix_path}: found {actual_bytes} bytes; expected {expected_bytes}"
        )

    samples = read_samples(sample_path, sample_id_columns)
    if any(len(values) != sample_count for values in samples.values()):
        raise ValueError(f"{sample_path}: sample count differs from score metadata")
    score_columns = read_score_columns(score_metadata_path, score_count)
    duplicate_names = set(sample_id_columns).intersection(score_columns)
    if duplicate_names:
        raise ValueError(
            "sample and score columns overlap: " + ", ".join(sorted(duplicate_names))
        )

    matrix = np.memmap(
        matrix_path,
        dtype=np.dtype("<f8"),
        mode="r",
        shape=(score_count, sample_count),
        order="C",
    )
    fields = [pa.field(name, pa.string(), nullable=False) for name in sample_id_columns]
    fields.extend(pa.field(name, pa.float64(), nullable=False) for name in score_columns)
    schema = pa.schema(fields).with_metadata(
        {
            b"pgensparsescore_format": EXPECTED_FORMAT.encode(),
            b"pgensparsescore_layout": b"sample-rows,named-score-columns",
        }
    )

    temporary_path = Path(str(output_path) + ".tmp")
    temporary_metadata_path = Path(str(output_path) + ".json.tmp")
    output_metadata_path = Path(str(output_path) + ".json")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path.unlink(missing_ok=True)
    temporary_metadata_path.unlink(missing_ok=True)
    writer: pq.ParquetWriter | None = None
    try:
        writer = pq.ParquetWriter(
            temporary_path,
            schema,
            compression=compression,
            compression_level=compression_level,
            use_dictionary=sample_id_columns,
            write_statistics=sample_id_columns,
        )
        for begin in range(0, sample_count, batch_rows):
            end = min(begin + batch_rows, sample_count)
            arrays = [
                pa.array(samples[name][begin:end], type=pa.string())
                for name in sample_id_columns
            ]
            arrays.extend(
                pa.array(matrix[index, begin:end], type=pa.float64())
                for index in range(score_count)
            )
            writer.write_batch(pa.RecordBatch.from_arrays(arrays, schema=schema))
        writer.close()
        writer = None

        result = {
            "format": "pgensparsescore-wide-parquet-v1",
            "path": output_path.name,
            "sample_id_columns": sample_id_columns,
            "sample_rows": sample_count,
            "score_columns": score_count,
            "compression": compression,
            "compression_level": compression_level,
            "row_group_rows": batch_rows,
        }
        with temporary_metadata_path.open("w", encoding="utf-8") as handle:
            json.dump(result, handle, indent=2)
            handle.write("\n")
        os.replace(temporary_path, output_path)
        os.replace(temporary_metadata_path, output_metadata_path)
        return result
    finally:
        if writer is not None:
            writer.close()
        temporary_path.unlink(missing_ok=True)
        temporary_metadata_path.unlink(missing_ok=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--score-metadata", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--batch-rows", type=int, default=4096)
    parser.add_argument("--compression", default="zstd")
    parser.add_argument("--compression-level", type=int, default=3)
    args = parser.parse_args()
    result = convert(
        args.metadata,
        args.score_metadata,
        args.output,
        args.batch_rows,
        args.compression,
        args.compression_level,
    )
    print(
        f"wrote {result['sample_rows']} sample rows and "
        f"{result['score_columns']} named score columns to {args.output}"
    )


if __name__ == "__main__":
    main()
