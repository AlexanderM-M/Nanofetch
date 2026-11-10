# Contributing

For bugs and feature requests, use the
[issue templates](https://github.com/AlexanderM-M/Nanofetch/issues/new/choose).
Check the [troubleshooting guide](docs/troubleshooting.md) first. Include the
version, command, environment, and expected behavior; use synthetic examples
instead of private alignment data.

## Development setup

From the repository root:

Create a virtual environment, install the editable test dependencies, and run the
test suite before opening a pull request:

```bash
python -m venv .venv
. .venv/bin/activate
python -m pip install -e '.[test]'
python -m nanofetch --version
python -m pytest
```

The editable install makes source changes available immediately. Run these
commands again when reusing an environment that cannot import `nanofetch`.

## Pull requests

Explain the problem, resulting behavior, and validation in the pull request.
Keep the README focused on getting started; detailed examples belong in the
[usage guide](docs/usage.md), and common errors in
[troubleshooting](docs/troubleshooting.md).

Annotation updates must be pinned, reproducible, and accompanied by updated
source URLs and SHA-256 digests in `annotations.json`. Do not edit the compressed
gene tables manually; regenerate them with `scripts/build_annotations.py`.

Please add tests for behavior changes. Clinical claims and opaque panel changes
are out of scope.

## Releases

Release artifacts must be built from a clean checkout or after removing ignored
`build/`, `dist/`, and `*.egg-info` directories. Install and smoke-test the built
wheel itself before uploading it; source-tree tests do not verify wheel contents.
