# pgensparsescore

`pgensparsescore` scores many polygenic scores in one pass over the relevant
variants of a PGEN file.  Its weight index is variant-major: a decoded dosage
vector is applied only to the scores which contain that variant.

This is an early correctness-first implementation.  It supports biallelic
variants, PGEN hardcalls and dosages, mean imputation of missing dosages, and
weights whose effect allele is either REF or ALT.  The allele order in a
weight row does not need to match the PVAR allele order.  PGEN dosages are used
on the stored 0--2 scale; chromosome- and sex-specific ploidy transformations
are not yet performed.

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
be the PVAR REF/ALT pair, but they may appear in either order.  For example,
both of these are valid for a PVAR record with `REF=C` and `ALT=T`:

```text
EFFECT_ALLELE=T  OTHER_ALLELE=C
EFFECT_ALLELE=C  OTHER_ALLELE=T
```

The scorer does not perform strand-complement guessing.  A matched variant
whose two alleles do not exactly equal the PVAR REF/ALT pair is an error.
Multiallelic PVAR records are rejected when referenced by a score.

## Run

```sh
pgensparsescore \
  --pgen cohort.pgen \
  --pvar cohort.pvar \
  --psam cohort.psam \
  --manifest scores.tsv \
  --out results/catalog
```

The primary output, `catalog.scores.tsv.gz`, is a gzip-compressed wide table.
Each row is a sample, followed by one column named for every manifest `SCORE`.
`IID` is required in the PSAM.  If the PSAM also has `FID`, it is preserved;
the scorer does not invent one when it is absent:

```text
IID\tscore_a\tscore_b
sample1\t0.125\t-1.75
sample2\t0.5\t0.25
```

Additional outputs are:

- `catalog.score-metadata.tsv`: score order and matching/QC counts.
- `catalog.json`: output dimensions and run-level counters.

The scorer uses a file-backed score-major matrix while applying variants,
since that is the efficient update layout.  It transposes that working matrix
in bounded-memory blocks when writing the sample-major table, then removes the
working file.  Score-major layout is therefore an implementation detail, not
part of the output contract.

## Scoring semantics

The input allele order determines which of two equivalent calculations is
used; it is not an input restriction:

| Weight row | Contribution |
| --- | --- |
| effect = PVAR ALT, other = PVAR REF | `w * ALT_DOSAGE` |
| effect = PVAR REF, other = PVAR ALT | `w * (2 - ALT_DOSAGE)` |

Internally, the second expression is rewritten as
`2 * w - w * ALT_DOSAGE`.  This permits the scoring kernel to operate entirely
on ALT dosages without changing the meaning of a REF-effect weight.  Missing
ALT dosages are mean-imputed from nonmissing samples at that variant.  An
all-missing scored variant is an error.

When `pgenlib` can return a common dosage plus a short difference list, the
common contribution is accumulated once per score and only the differing
samples are updated.  Otherwise, the dense dosage vector is decoded and
applied.  The score matrix is file-backed, so RAM does not grow as
`sample_count * score_count`.

## License

The application is GPL-3.0-only.  `pgenlib` is a separate LGPL-3.0 dependency.
