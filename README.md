# pgensparsescore

`pgensparsescore` scores many polygenic scores in one pass over the relevant
variants of a PGEN file.  Its weight index is variant-major: a decoded dosage
vector is applied only to the scores which contain that variant.

This is an early correctness-first implementation.  It supports biallelic
variants, PGEN hardcalls and dosages, fixed-frequency or cohort-frequency mean
imputation, and weights whose effect allele is either REF or ALT.  The allele
order in a weight row does not need to match the PVAR allele order.  PGEN
dosages are used on the stored 0--2 scale; chromosome- and sex-specific ploidy
transformations are not yet performed.

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
Development packages for zlib and zstd are also required.

## Inputs

The score manifest is a tab-separated file with `SCORE` and `PATH` columns.
Relative paths are resolved relative to the manifest.

```text
SCORE\tPATH
score_a\tweights/score_a.tsv.gz
score_b\tweights/score_b.tsv.gz
```

Each weight file is a tab-separated plain-text, gzip-compressed, or
zstd-compressed file with these columns:

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

To score variant-sharded PGENs, such as one file per chromosome, provide a
tab-separated PGEN list instead:

```text
PGEN\tPVAR\tPSAM
cohort.chr1.pgen\tcohort.chr1.pvar.zst\tcohort.chr1.psam
cohort.chr2.pgen\tcohort.chr2.pvar.zst\tcohort.chr2.psam
```

```sh
pgensparsescore \
  --pfile-list cohort.pfiles.tsv \
  --manifest scores.tsv \
  --out results/catalog
```

Relative paths in the PGEN list are resolved relative to that list.  Every
PGEN must represent a distinct set of variants and have exactly the same
samples, FIDs when present, and sample order.  The score manifest is compiled
once against the union of the PVARs, then partitioned by PGEN.  PGENs are
decoded sequentially into one file-backed accumulator and the output table is
written only once.

### Fixed-frequency imputation

Use `--read-freq` when the same reference frequency must impute missing
genotypes in more than one target dataset:

```sh
pgensparsescore \
  --pfile-list cohort.pfiles.tsv \
  --manifest scores.tsv \
  --read-freq reference.acount.zst \
  --error-on-missing-freq \
  --out results/catalog
```

The frequency file can be plain text, gzip, or zstd.  It must have `ID`,
`REF`, and `ALT` (or `ALT1`), plus one of the following representations:

- `ALT_DOSAGE_MEAN`;
- `ALT_FREQS` or `ALT1_FREQ`;
- `ALT_CTS` or `ALT1_CT`, together with `OBS_CT`; or
- `REF_CT`, together with `OBS_CT`.

This includes biallelic PLINK 2 `.acount[.zst]` and `.afreq[.zst]` files.  An
allele frequency `p` is converted to expected ALT dosage `2p`.  For every
scored variant found in the frequency table, the frequency REF and ALT must
exactly match the PVAR REF and ALT; disagreement is always an error.

With `--error-on-missing-freq`, every scored variant must have a frequency
row.  Without that flag, variants absent from the frequency table fall back to
the nonmissing mean in the PGEN currently being scored, and the fallback is
counted in the JSON metadata.  For cross-dataset comparability, use the strict
form shown above.

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
ALT dosages are mean-imputed from the external expected ALT dosage when
`--read-freq` supplies the variant.  Otherwise they are imputed from nonmissing
samples at that variant.  An all-missing scored variant is permitted with an
external frequency and is an error without one.

When `pgenlib` can return a common dosage plus a short difference list, the
common contribution is accumulated once per score and only the differing
samples are updated.  Otherwise, the dense dosage vector is decoded and
applied.  The score matrix is file-backed, so RAM does not grow as
`sample_count * score_count`.

## License

The application is GPL-3.0-only.  `pgenlib` is a separate LGPL-3.0 dependency.
