# PyFEM documentation

Sphinx sources in this directory. Build after installing the package:

```bash
uv sync --extra docs --no-dev
uv run sphinx-build -M html doc doc/_build
```

Output: `doc/_build/html/index.html`

`--no-dev` matches CI and Read the Docs. For a single environment with pytest/ruff as well, use `uv sync --extra docs` without `--no-dev`.

API reference is generated at build time by `sphinx-autoapi` under `doc/api/` (gitignored; not committed).

## Hosting

- **Read the Docs (primary):** https://pyfem.readthedocs.io/
- **GitHub Pages (mirror):** built on every PR and push; published only on **push to `main`** via `.github/workflows/ci.yml`. Requires repository **Settings → Pages → Build and deployment → Source: GitHub Actions**.

## Linux note

Doc builds may require `libgl1` because Sphinx viewcode imports VTK-related modules. CI and RTD install it via apt.

Optional live reload (no extra project dependency):

```bash
uv run --with sphinx-autobuild sphinx-autobuild doc doc/_build/html --watch pyfem --re-ignore "_build/*"
```
