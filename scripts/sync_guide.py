"""Build the user-guide pages from the exuber vignettes.

The vignettes are the single source. Each guide page is one or more vignettes,
knitted with real output, stripped of their YAML header, and written to
src/content/guide/<slug>.md with figures in public/guide-figs/. Requires R with
exuber and knitr installed (set R_LIBS if exuber lives in a private library).

    python scripts/sync_guide.py
"""
import re
import subprocess
import sys
import tempfile
from pathlib import Path

WEB = Path(__file__).resolve().parent.parent
VIGNETTES = WEB.parent / "exuber" / "vignettes"
OUT = WEB / "src" / "content" / "guide"
FIGS = WEB / "public" / "guide-figs"

# slug, title, blurb, group, order, replication family (or None), vignettes
PAGES = [
    ("pipeline", "Results, tidying and plotting", "How a test result flows through summary(), diagnostics(), datestamp(), tidy() and autoplot(), and what the function names mean.", "Guide", 4, None, ["naming-and-analysis", "plotting"]),
    ("volatility-robustness", "Volatility-robustness tests", "What to do when the variance of your series changes over time and the standard critical values over-reject.", "Beyond PSY", 1, "volatility-robustness", ["volatility-robust-radf", "radf-tt"]),
    ("dating", "Dating and root inference", "Date a bubble you already believe is there, and estimate how fast it grows.", "Beyond PSY", 2, "dating-and-root-inference", ["dating-methods", "root-inference", "experimental-methods"]),
    ("monitoring", "Real-time monitoring for bubbles", "Watch new observations arrive and raise an alarm the first time a boundary is crossed.", "Beyond PSY", 3, "monitoring", ["monitoring"]),
    ("multivariate", "Multivariate bubble tests", "Test whether two series bubble together, and whether a bubble in one spills into another.", "Beyond PSY", 4, "multivariate", ["co-explosivity"]),
    ("alternatives", "Alternative paradigms", "Tests built on different ideas from the recursive ADF: locally best invariant, stochastic-coefficient and quantile tests.", "Beyond PSY", 5, "alternative-paradigms", ["alternative-tests"]),
    ("simulation", "Simulating bubbles", "Generate bubble data to study size, power and detection delay.", "Beyond PSY", 6, "simulation-dgps", ["simulation"]),
]

SLUG_OF = {v: p[0] for p in PAGES for v in p[6]}
SLUG_OF["exuber"] = None  # the getting-started page


def knit(vignette: str, tmp: Path) -> str:
    md = tmp / f"{vignette}.md"
    script = (
        "library(knitr);"
        f"opts_chunk$set(fig.path='{tmp}/figs/{vignette}-',dev='svg',fig.width=7,fig.height=3.6);"
        f"knit('{VIGNETTES / (vignette + '.Rmd')}',output='{md}',quiet=TRUE)"
    )
    subprocess.run(["Rscript", "-e", script], check=True, cwd=tmp)
    return md.read_text(encoding="utf-8")


def clean(text: str, vignette: str) -> tuple[str, str]:
    m = re.match(r"---\n(.*?)\n---\n", text, re.S)
    title = re.search(r'title:\s*"(.*?)"', m.group(1)).group(1) if m else vignette
    body = text[m.end():] if m else text
    body = re.sub(r"\n{3,}", "\n\n", body).strip()
    # vignette("x") -> link to the matching guide page
    def link(mo):
        slug = SLUG_OF.get(mo.group(1))
        if slug:
            return f"[the {PAGE_TITLE[slug].lower()} page](/guide/{slug})"
        return "[Getting started](/guide)" if mo.group(1) == "exuber" else mo.group(0)
    body = re.sub(r"`?vignette\(\"([\w-]+)\"\)`?", link, body)
    body = re.sub(r"!\[([^\]]*)\]\(.*?figs/([^)]+)\)", r"![\1](/guide-figs/\2)", body)
    return title, body


PAGE_TITLE = {p[0]: p[1] for p in PAGES}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as t:
        tmp = Path(t)
        (tmp / "figs").mkdir()
        for slug, title, blurb, group, order, family, vigs in PAGES:
            parts = []
            for v in vigs:
                vt, body = clean(knit(v, tmp), v)
                if len(vigs) > 1:
                    body = re.sub(r"^(#{2,5}) ", lambda m: "#" * (len(m.group(1)) + 1) + " ", body, flags=re.M)
                    body = f"## {vt}\n\n{body}"
                parts.append(body)
            front = (
                f'---\ntitle: "{title}"\nblurb: "{blurb}"\ngroup: "{group}"\n'
                f'order: {order}\nfamily: "{family or ""}"\n---\n'
            )
            (OUT / f"{slug}.md").write_text(front + "\n\n".join(parts) + "\n", encoding="utf-8")
            print(f"guide   : {slug} <- {', '.join(vigs)}")
        for f in (tmp / "figs").glob("*.svg"):
            (FIGS / f.name).write_bytes(f.read_bytes())
        for f in (tmp / "figs").glob("*.png"):
            (FIGS / f.name).write_bytes(f.read_bytes())


if __name__ == "__main__":
    sys.exit(main())
