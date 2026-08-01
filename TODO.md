# TODO

## Current priorities

1. **Add a deterministic secondary tie-break to candidate-region peak ranking**
   - Where: `ABC/workflow/scripts/peaks.py` lines 37 and 83, both currently `"sort -nr -k 4"` (descending by accessibility count, column 4), immediately followed by `head -n {nStrongestPeaks}` to take the top-N peaks.
   - Problem: when many peaks share the exact same count at the cutoff boundary (observed: tie groups of hundreds to thousands of peaks in real full-genome data), GNU `sort`'s tie-break order for equal keys is unspecified and can vary across environments/invocations, so which peaks land just inside vs. just outside the top-N is nondeterministic — even with byte-identical input. Root-caused and confirmed not a container/conda bug; see `TESTING.md`'s "Known test-coverage gaps" for the summary and git history (`dev` branch, commit `88f5dfa` and earlier) for the full investigation writeup.
   - Suggested fix: extend the sort key so no two rows can tie, e.g. `sort -k4,4nr -k1,1 -k2,2n -k3,3n` (count descending, then chr/start/end ascending as a fully-resolving secondary key on genomic coordinates). Since no two peaks share the same chr/start/end, this makes the sort's result fully deterministic regardless of multi-threading or `sort` implementation/version — no need for `--parallel=1` or `LC_ALL` pinning once ties are fully broken by the key itself.
   - Verify: rerun a chr22 and/or full-genome conda-vs-singularity comparison after the fix and confirm `candidateRegions.bed` is now byte-identical run-to-run (not just same total count).
   - Low urgency: effect on final threshold-filtered predictions is negligible (single-digit pair differences in the confident set observed), so this is about reproducibility/determinism, not correctness.
   - **Blocked on an external dependency**: `ABC/` is a nested submodule pointing at `broadinstitute/ABC-Enhancer-Gene-Prediction.git` (not an EngreitzLab repo), itself nested inside the `ENCODE_rE2G` submodule. Landing this fix durably requires: (1) fork/patch upstream ABC (or get a PR merged there), (2) bump `ENCODE_rE2G`'s ABC submodule pointer to the fixed commit, (3) bump scE2G's `ENCODE_rE2G` submodule pointer to that updated commit. A local-only edit can verify the fix works but can't be committed/shared normally, since we don't own that upstream repo.

2. **Review code documentation, organization, and clarity**
   - Audit inline comments and docstrings
   - Review README and user-facing documentation
   - Assess script organization and naming conventions
   - Identify areas needing clarification

3. **Merge `confidence-intervals-cell-type-comparison` work into `dev`/`main`**
   - CI subsampling workflow (`workflow/Snakefile_ci`, `workflow/rules/ci_*.smk`, `workflow/scripts/ci/`) and the cross-cell-type comparison/visualization scripts (`workflow/scripts/visualization/`) are already implemented — see `CLAUDE.md` for the full description — but live only on `confidence-intervals-cell-type-comparison` (which itself depends on `flatten-encode-re2g-submodule`'s flattened submodule namespace).
   - Remaining work: merge `flatten-encode-re2g-submodule` forward first (submodule flattening + checkpoint removal), then rebase/merge the CI-comparison branch on top, then land both on `dev`/`main`.

## Test data

WTC11 differentiation timecourse (iPSC → EC, 5 cell types): `results/2026_0210_WTC11_EC/` — see `CLAUDE.md` for details. Precomputed comparison outputs already exist there for exercising the cross-cell-type comparison scripts once item 3 lands.
