# Contributing

Thanks for your interest in improving this package. This guide covers how to report
problems, set up a development environment, run the tests, and propose changes.

By participating, you agree to abide by the [Code of Conduct](CODE_OF_CONDUCT.md).

## Reporting issues

Open a [GitHub issue](https://github.com/gynecoloji/Annotation_st/issues) with:

- the `annot-st` command (or Python call) you ran and the config file, with marker-file
  paths if any,
- what you expected vs. what happened,
- the stage log from `<outdir>/logs/`, `annot-st status -c <config>` output, and the
  versions printed by `pip show annotation_st scanpy squidpy anndata`.

## Development setup

```bash
git clone https://github.com/gynecoloji/Annotation_st.git
cd annotation_st
python -m venv .venv && source .venv/bin/activate     # or a conda env with python >= 3.10
pip install -e ".[test]"
pip install spatioloji-s                               # optional: the spatioloji_s bridge tests
```

On easley the `spatioloji` conda environment already has every dependency:
`/users/jiwang1/.conda/envs/spatioloji/bin/python -m pip install -e . --no-deps --no-build-isolation`.

## Running the tests

```bash
MPLBACKEND=Agg python -m pytest tests -q          # about 3 minutes; synthetic data, no network
```

The tests build a synthetic Xenium-like dataset with planted populations (`tests/conftest.py`)
and check every decision rule, the colour scheme, staleness tracking, the command line and
the spatioloji_s bridge. Tests that need `spatioloji_s` skip themselves when it is not installed.

The default marker library is data (`src/annotation_st/data/marker_library/`). Changing it
means re-running `tools/build_marker_references.py` (needs network access): every gene must
keep a reference, and a test enforces that.

## Proposing changes

1. Branch off `main`; keep one topic per pull request.
2. Write the failing test first, then the change; keep the suite green.
3. Use [Conventional Commits](https://www.conventionalcommits.org/) in commit messages
   (`feat:`, `fix:`, `docs:`, `refactor:`, `test:`, `chore:`). release-please turns them into
   the changelog and the next version number.
4. If a change alters labels on the validated OV section, say so in the PR and update
   `validation/README.md` after re-running `validation/run_ov_validation.sbatch`.
