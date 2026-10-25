"""Build the Kaggle notebook from the percent-format sources in notebook/src/*.py.

Each source file is plain Python with `# %%` cell markers (`# %% [markdown]` for text cells),
so every code cell can be syntax-checked before the notebook is pushed.

    python notebook/build_nb.py            # writes notebook/kaggle/<slug>.ipynb
"""
import ast
import json
import sys
from pathlib import Path

import nbformat

HERE = Path(__file__).resolve().parent
SRC = HERE / "src"
OUT_DIR = HERE / "kaggle"


def split_cells(text):
    cells, kind, buf = [], None, []
    for line in text.splitlines():
        if line.startswith("# %%"):
            if kind is not None:
                cells.append((kind, buf))
            kind = "markdown" if "[markdown]" in line else "code"
            buf = []
        else:
            if kind is None:        # preamble before the first marker is ignored
                continue
            buf.append(line)
    if kind is not None:
        cells.append((kind, buf))
    return cells


def main():
    meta = json.loads((OUT_DIR / "kernel-metadata.json").read_text())
    nb = nbformat.v4.new_notebook()
    nb.metadata["kernelspec"] = {"display_name": "Python 3", "language": "python", "name": "python3"}
    nb.metadata["language_info"] = {"name": "python"}
    n_code = n_md = 0
    for src in sorted(SRC.glob("*.py")):
        for kind, buf in split_cells(src.read_text(encoding="utf-8")):
            while buf and not buf[0].strip():
                buf.pop(0)
            while buf and not buf[-1].strip():
                buf.pop()
            if not buf:
                continue
            if kind == "markdown":
                body = "\n".join(l[2:] if l.startswith("# ") else l.lstrip("#") for l in buf)
                nb.cells.append(nbformat.v4.new_markdown_cell(body))
                n_md += 1
            else:
                body = "\n".join(buf)
                try:
                    ast.parse(body)
                except SyntaxError as e:
                    sys.exit(f"SyntaxError in {src.name}: line {e.lineno}: {e.msg}\n{e.text}")
                nb.cells.append(nbformat.v4.new_code_cell(body))
                n_code += 1
    nbformat.validate(nb)
    out = OUT_DIR / meta["code_file"]
    nbformat.write(nb, out)
    print(f"wrote {out} ({n_code} code cells, {n_md} markdown cells, {out.stat().st_size / 1024:.0f} KB)")


if __name__ == "__main__":
    main()
