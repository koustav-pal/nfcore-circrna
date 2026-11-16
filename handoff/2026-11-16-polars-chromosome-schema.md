# Polars chromosome schema fix — 2026-11-16

## Scope and state

Repository: `koustav-pal/nfcore-circrna`, branch `dev`.
Request: implement `prompts/fix_polars_schema_bug.md`, push changes, save handoff here.
Starting commit: `b279622f`. Existing untracked `prompts/` preserved and excluded from commit.
No repository or ancestor AGENTS.md/CLAUDE.md instructions found.

Implementation and regression coverage are included in the commit containing this handoff.
See that commit's diff and `CHANGELOG.md` for the file inventory; do not reimplement.

## Diagnosis and decisions

Polars infers integer chromosomes from early numeric rows and then fails on late `X`.
Explicit chromosome string schemas fix the observed parsing failure while retaining numeric coordinates.
All 11 affected Polars calls across combinebeds, jccirc, CIRI build_list, and circtools now override only chromosome dtype.
CIRI consensus uses `column_0`, confirmed against pinned Polars 2.0.0's headerless naming;
its existing projection, renaming and BED-to-annotation coordinate conversion remain unchanged.
CIRI-vis explicitly reads `Chr` as pandas string, including numeric-only input files.

Repository-wide reader audit found no additional affected reads:
- bed2gtf, CIRI3 coordinate readers, and circtools supporting-read schema already type chromosomes.
- majority_vote, readlength, matrix/join_samples, and PSIRC inputs have no separate chromosome column.
- matrix/join_samples call sites use identifier/gene metadata, not coordinate columns.

## Validation and continuation

Temporary environment: `/private/tmp/circrna-schema-venv` (not part of repository).
Pinned validation dependencies: Polars 2.0.0, pandas 2.3.3.

```sh
/private/tmp/circrna-schema-venv/bin/python -m unittest discover -s tests/python -v
git diff --check
```

Both unittest methods pass: 24 reader cases (12 readers × numeric-only/mixed input), plus
full CIRI build_list execution at max_shift 0 and 1 with BED12 consensus and annotation joins.
Assertions preserve chromosome strings, late `X`, integer start/end columns and expected coordinates.
The initial reader regression failed before the fix with chromosome dtype mismatches and late-X parse errors.
An initial reader-only run also passed on Polars 2.1.0 / pandas 3.1.1.

Native module regression added to `modules/local/combinebeds/filter/tests/main.nf.test`,
tagged `chromosome_regression`, following existing bed2gtf fixture conventions.
Nextflow and nf-test executables are unavailable locally; native nf-test and full pipeline execution remain unverified.
On a configured runner:

```sh
nf-test test modules/local/combinebeds/filter/tests/main.nf.test --tag chromosome_regression --profile docker
```

A direct full-filter before/after harness is available temporarily at
`/private/tmp/verify_circrna_filter.py`; it renders the real template, supplies 101 numeric
chromosomes followed by `X`. This passed: baseline failed with ``could not parse `X` as dtype `i64` ``;
fixed template emitted all 102 chromosomes with coordinates unchanged.

## Suggested skills

- `caveman`: user requested terse communication; remains active.
- `diagnose`: if subsequent pipeline failures arise, reproduce the exact failing reader first.
- `handoff`: update this document if a later session changes implementation or validation state.

No further code changes planned. Run native regression/CI and resume the user's pipeline to validate real data.
