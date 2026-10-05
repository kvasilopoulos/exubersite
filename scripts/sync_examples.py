"""Run every /reference example and keep its printed output and plots.

The code comes from src/data/reference.json (so run sync_content.py first).
Each topic is knitted against the working tree of exuber/ (pkgload::load_all),
with output interleaved after the code the way pkgdown shows it. Results go to
src/data/reference-output.json and figures to public/reference-figs/.

    python scripts/sync_examples.py [topic ...]
"""
import json
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

WEB = Path(__file__).resolve().parent.parent
PKG = WEB.parent / "exuber"
DATA = WEB / "src" / "data"
FIGS = WEB / "public" / "reference-figs"
RSCRIPT = os.environ.get("RSCRIPT", "Rscript")  # Windows: path to Rscript.exe

R = """
args <- commandArgs(TRUE)
suppressMessages(pkgload::load_all('%s', quiet = TRUE))
library(knitr)
knitr::opts_chunk$set(comment = '#>', collapse = TRUE, error = TRUE, dev = 'svg',
                      fig.width = 7, fig.height = 3.6)
knitr::opts_chunk$set(fig.path = paste0(args[2], '/'))
knitr::knit(args[1], output = sub('[.]Rmd$', '.md', args[1]), quiet = TRUE)
""" % PKG.as_posix()


def blocks(md: str, topic: str) -> list[dict]:
    """Split knitr markdown into code (with #> output) and image blocks."""
    out = []
    for part in re.split(r"(!\[[^\]]*\]\([^)]*\))", md):
        img = re.match(r"!\[[^\]]*\]\(.*?([^/)]+)\)", part)
        if img:
            out.append({"img": f"/reference-figs/{topic}/{img.group(1)}"})
            continue
        code = "\n".join(re.findall(r"``` ?r?\n(.*?)\n```", part, re.S)).strip()
        if code:
            out.append({"code": code})
    return out


def main() -> None:
    topics = json.loads((DATA / "reference.json").read_text(encoding="utf-8"))["topics"]
    want = sys.argv[1:] or [t for t, d in topics.items() if d.get("examples")]
    result = {}
    with tempfile.TemporaryDirectory() as t:
        tmp = Path(t)
        (tmp / "run.R").write_text(R, encoding="utf-8")
        for topic in want:
            rmd = tmp / f"{topic}.Rmd"
            rmd.write_text("```{r}\n" + topics[topic]["examples"] + "\n```\n", encoding="utf-8")
            figdir = tmp / "figs" / topic
            figdir.mkdir(parents=True, exist_ok=True)
            p = subprocess.run([RSCRIPT, str(tmp / "run.R"), str(rmd), figdir.as_posix()],
                               cwd=tmp, capture_output=True, text=True)
            md = rmd.with_suffix(".md")
            if p.returncode or not md.exists():
                print(f"examples: {topic} FAILED\n{p.stderr[-400:]}")
                continue
            result[topic] = blocks(md.read_text(encoding="utf-8"), topic)
            dest = FIGS / topic
            dest.mkdir(parents=True, exist_ok=True)
            for f in figdir.glob("*.svg"):
                (dest / f.name).write_bytes(f.read_bytes())
            print(f"examples: {topic} ({sum('img' in b for b in result[topic])} figures)")
    path = DATA / "reference-output.json"
    old = json.loads(path.read_text(encoding="utf-8")) if path.exists() and sys.argv[1:] else {}
    path.write_text(json.dumps({**old, **result}, indent=1), encoding="utf-8")


if __name__ == "__main__":
    sys.exit(main())
