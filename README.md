# pgensparsescore

`pgensparsescore` scores many polygenic scores in one pass over the relevant
variants of a PGEN file.  Its weight index is variant-major: a decoded dosage
vector is applied only to the scores which contain that variant.

This is an early correctness-first implementation.  It supports biallelic
variants, PGEN hardcalls and dosages, mean imputation of missing dosages, and
REF/ALT effect-allele orientation.  PGEN dosages are used on the stored
0--2 scale; chromosome- and sex-specific ploidy transformations are not yet
performed.

## Build

The executable links directly to the LGPL-3.0 `pgenlib` sources.  CMake can
fetch a pinned public revision of [plink-ng](https://github.com/chrchang/plink-ng),
or an existing checkout can be supplied for offline builds.

```sh
cmake -S . -B build -DPGENLIB_SOURCE_DIR=/path/to/plink-ng/2.0/include
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Omit `PGENLIB_SOURCE_DIR` to let CMake fetch the pinned revision.

## Inputs

The score manifest is a tab-separated file with `SCORE` and `PATH` columns.
Relative paths are resolved relative to the manifest.

```text
SCORE\tPATH
score_a\tweights/score_a.tsv.gz
score_b\tweights/score_b.tsv.gz
```

Each weight file is a tab-separated plain-text or gzip-compressed file with
these columns:

```text
SNP\tEFFECT_ALLELE\tOTHER_ALLELE\tEFFECT_ALLELE_WEIGHT
```

`SNP` is matched exactly to the PVAR `ID`.  The effect and other alleles must
match its REF/ALT pair in either order.  Multiallelic PVAR records are rejected
when referenced by a score.

## Run

```sh
pgensparsescore \
  --pgen cohort.pgen \
  --pvar cohort.pvar \
  --psam cohort.psam \
  --manifest scores.tsv \
  --out results/catalog
```

The output is deliberately simple and language-independent:

- `catalog.scores.bin`: little-endian float64 matrix in C order, with shape
  `(score_count, sample_count)`.
- `catalog.scores.tsv`: score order and matching/QC counts.
- `catalog.samples.tsv`: sample order.
- `catalog.json`: matrix shape, dtype, order, and run-level counters.

For example, Python can open the matrix without copying it:

```python
scores = numpy.memmap(
    "results/catalog.scores.bin",
    dtype="<f8",
    mode="r",
    shape=(score_count, sample_count),
)
```

## Scoring semantics

For a weight `w`, an ALT effect allele contributes `w * ALT_DOSAGE`.  A REF
effect allele contributes `w * (2 - ALT_DOSAGE)`.  Missing ALT dosages are
mean-imputed from nonmissing samples at that variant.  An all-missing scored
variant is an error.

When `pgenlib` can return a common dosage plus a short difference list, the
common contribution is accumulated once per score and only the differing
samples are updated.  Otherwise, the dense dosage vector is decoded and
applied.  The score matrix is file-backed, so RAM does not grow as
`sample_count * score_count`.

## License

The application is GPL-3.0-only.  `pgenlib` is a separate LGPL-3.0 dependency.
