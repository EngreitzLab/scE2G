# Handover

Purpose: let a new Claude Code session pick up exactly where the last one left off if this
session is killed/restarted. Read this alongside `CLAUDE.md` (project overview, branch map)
and `TODO.md` (priority list). Keep this file updated as work progresses — see instruction
in `CLAUDE.md`.

## Current branch: `dev`

**Update 2026-08-01: repo cleanup after PR #109 (`add-chr22-and-unit-tests`) merged to `main`.**

- Deleted 9 local + 2 remote branches whose commits were already fully merged into `main`
  (`add_stats`, `containerize`, `hover_plots`, `list`, `match_gene_names`, `maya_updates`,
  `megamap_v3_models`, `qc_plots`, `fix_gene_mapping`, `origin/613weilin-patch-1`,
  `origin/ms_conda_env`).
- Deleted `ci/step8a-path-filtering` (identical to old `dev` tip, no unique work) and
  `ci/testing-changes-to-main` (a prior attempt to merge `dev`'s full testing/CI overhaul into
  `main`, superseded by the smaller squashed `add-chr22-and-unit-tests` PR that actually landed).
- **`dev`'s testing setup was reset to match `main`'s.** `dev` had diverged from `main` before
  containerization (`#97`), `igv_rule_precedence` (`#101`), the jq/getcwd fixes (`#106`/`#107`),
  and had independently built a more elaborate tier-1/tier-2/tier-3 test suite plus CircleCI
  path-filtering, env-build caching, and nightly-rebuild plumbing — none of which ever landed on
  `main`. `main` instead got a simpler, squashed version (unit tests + chr22-scoped E2E test,
  single static CircleCI job, see `TESTING.md`) via `add-chr22-and-unit-tests`. Per user decision,
  `main`'s simpler setup is authoritative going forward: `dev` was rebuilt from `main`'s current
  tip, carrying forward only non-testing-related content (`CLAUDE.md`, `TODO.md`, `HANDOVER.md`,
  `.gitignore`'s `scratch/` entry). The old `dev`'s elaborate testing-infra commits and detailed
  investigation history (container-vs-conda bug hunt, peak tie-break root-cause) remain
  accessible via `git log` on the old `dev` tip (`fb91db9`) if ever needed for reference — not
  carried forward into this file since `TESTING.md` on `main` already documents the current
  conclusions.
- Branches intentionally left untouched (real in-progress work, not "random"):
  `confidence-intervals-cell-type-comparison`, `flatten-encode-re2g-submodule` (see `TODO.md`
  item 3), plus other collaborators' branches (`kendall_optimization`, `kaybrand-declare-nightly`,
  `mayasheth-patch-1`, `mouse_gene_annotation`, `rosa_mod`, `sparse_input`,
  `sparse_input_rebase`) — not this session's to delete.

## Not yet committed / pushed

This rebuild of `dev` is staged locally on branch `dev-new` (based on `origin/main` tip
`06f2792`), not yet force-pushed to `origin/dev`. **Standing instruction: ask before every `git
commit`/`git push`, and especially before a force-push that rewrites a shared branch's history.**

## Next steps

- Confirm the rebuilt `dev` looks right, then force-push `dev-new` over `origin/dev` (and update
  the local `dev` ref to match).
- Pick up `TODO.md` priorities — item 1 (peak tie-break determinism, blocked on upstream ABC) and
  item 3 (merging the confidence-interval/cross-cell-type-comparison work forward) are the two
  substantive open threads.
