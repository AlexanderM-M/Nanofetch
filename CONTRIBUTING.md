# Contributing

Create a virtual environment, install the editable test dependencies, and run the
test suite before opening a pull request:

```bash
python -m venv .venv
. .venv/bin/activate
python -m pip install -e '.[test]'
pytest
```

Annotation updates must be pinned, reproducible, and accompanied by updated
source URLs and SHA-256 digests in `annotations.json`. Do not edit the compressed
gene tables manually; regenerate them with `scripts/build_annotations.py`.

Please add tests for behavior changes. Clinical claims and opaque panel changes
are out of scope.

Release artifacts must be built from a clean checkout or after removing ignored
`build/`, `dist/`, and `*.egg-info` directories. Install and smoke-test the built
wheel itself before uploading it; source-tree tests do not verify wheel contents.
