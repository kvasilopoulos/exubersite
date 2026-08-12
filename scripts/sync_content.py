#!/usr/bin/env python
"""Vendors curated content from the project root into website/.

website/ is its own git repo and Netlify builds it from GitHub, but
docs/enhancements/ and exuber/ live at the project root, outside it --
nothing there exists at Netlify build time. So this runs manually and its
output is committed. No submodule, no CI step, no Netlify build change.

    python website/scripts/sync_content.py

Three jobs:
  1. curated docs/enhancements/*.md  -> src/content/replication/*.md
  2. docs/enhancements/replication/**/*.R -> src/replication-code/*.R
  3. exuber/_pkgdown.yml + man/*.Rd  -> src/data/reference.json,
     including a full parse of every Rd file (usage, arguments, value,
     examples, seealso, ...) so /reference/[topic] pages render natively
     instead of linking out to pkgdown.

Requires PyYAML (``pip install pyyaml``); everything else is stdlib.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import yaml

# ---- paths ------------------------------------------------------------------

WEB = Path(__file__).resolve().parent.parent
ROOT = WEB.parent
ENH = ROOT / "docs" / "enhancements"
EXUBER = ROOT / "exuber"

assert ENH.is_dir(), ENH
assert (EXUBER / "man").is_dir(), EXUBER / "man"

# Curated set. Deliberately withheld: papers/ (52 copyrighted PDFs),
# status.md (a conversational dump), README.md / SUMMARY.md (internal bundle
# log), open-research-directions.md and practitioner-guidance.md (backlog,
# not replication).
PAGES = [
    ("volatility-robustness", "Tests robust to time-varying innovation variance: time-transformed, kernel-purged, WLS, sign-based, and stochastic-coefficient routes."),
    ("dating-and-root-inference", "Origination, collapse and recovery dates, plus confidence intervals on the explosive root itself."),
    ("monitoring", "Sequential and real-time detection: training-vs-monitoring orchestration, CUSUM families, and closed-form boundaries."),
    ("multivariate", "Panel and cross-series tests -- common bubbles, co-bubbles, and bubble contagion."),
    ("alternative-paradigms", "Non-ADF-family approaches, principally the quantile-based global test and its recursive monitoring extension."),
    ("simulation-dgps", "Data-generating processes for the axes the original sim_*() functions do not cover."),
    ("references", "Full bibliography behind the replication record, organised by methodological family."),
]
PUBLISHED = {slug for slug, _ in PAGES}

# ============================================================================
# 1. markdown link rewriting (docs/enhancements/*.md -> site content)
# ============================================================================

LINK_RE = re.compile(r"\[[^\]]*\]\([^)\s]+\)")


def rewrite_one(m: re.Match) -> str:
    lnk = m.group(0)
    text = re.sub(r"^\[((?:.|\n)*)\]\([^)]*\)$", r"\1", lnk)
    target = re.sub(r"^\[(?:.|\n)*\]\(([^)]*)\)$", r"\1", lnk)

    if re.match(r"^(https?:|mailto:|#)", target):
        return lnk

    bare = re.sub(r"^(\.\./)+", "", target)
    anchor = re.search(r"#.*$", bare)
    anchor = anchor.group(0) if anchor else ""
    path = re.sub(r"#.*$", "", bare)

    if path == "":
        return lnk
    if path.startswith("papers/"):
        return text
    if path.endswith(".md"):
        slug = Path(path).stem
        if slug in PUBLISHED:
            return f"[{text}](/replication/{slug}{anchor})"
        return text
    if re.match(r"^replication/.*\.R$", path):
        return f"[{text}](#script-{Path(path).stem})"
    return text


def rewrite_links(txt: str) -> str:
    return LINK_RE.sub(rewrite_one, txt)


# Bare prose references to withheld or project-root paths -- "Local copy:
# papers/x/y.pdf.", "open (institutional access) -- papers/x/y.pdf",
# "`docs/enhancements/replication/fam/script.R`". Internal bookkeeping about
# where a local file sits; means nothing to a visitor and points at files the
# site does not ship.
_LOCAL_PATH_RE = re.compile(
    r"(?:[,;]?\s*(?:[Ll]ocal\s+cop(?:y|ies):?|[Ll]ocal\s+access:?|[Cc]opy:?|[Aa]ccess:?)?\s*(?:[-\u2013\u2014]\s*)?)"
    r"papers/\S+?\.(?:pdf|txt)[,.]?"
)


def strip_local_paths(txt: str) -> str:
    txt = _LOCAL_PATH_RE.sub("", txt)
    txt = re.sub(r"docs/enhancements/replication/[a-z-]+/([A-Za-z0-9_]+)\.R", r"\1.R", txt)
    txt = re.sub(r"docs/enhancements/([a-z-]+)\.md", r"/replication/\1", txt)
    txt = re.sub(r"\s*PDFs: `papers/`\.", "", txt)
    txt = txt.replace("`papers/`", "the paper library")
    return re.sub(r"[ \t]+\n", "\n", txt)


# references.md titles its sections "## Multivariate -- [multivariate.md](...)".
# The trailing link pushes the heading slug to `multivariate--multivariatemd`,
# which is not what the sibling files link to. Drop it so the slug is just the
# section name; scripts/check-links.mjs verifies every anchor after a build.
_HEADING_LINK_RE = re.compile(
    r"(?m)^(#{2,3} [^\n\[]+?)\s+[\u2014\u2013-]\s+(?:\[[^\]]+\]\([^)]*\)|[a-z0-9-]+\.md)\s*$"
)


def normalise_headings(txt: str) -> str:
    txt = _HEADING_LINK_RE.sub(r"\1", txt)
    # Stale in the source: alternative-paradigms.md still points at the
    # anchor references.md used before it was reorganised out of tier order.
    return txt.replace("#tier-7--quantile-based-detection", "#alternative-paradigms")


def yaml_str(x: str) -> str:
    return '"' + x.replace("\\", "\\\\").replace('"', '\\"') + '"'


def sync_markdown() -> None:
    out_dir = WEB / "src" / "content" / "replication"
    if out_dir.exists():
        for f in out_dir.glob("*.md"):
            f.unlink()
    out_dir.mkdir(parents=True, exist_ok=True)

    for order, (slug, blurb) in enumerate(PAGES, start=1):
        src = ENH / f"{slug}.md"
        assert src.exists(), src
        lines = src.read_text(encoding="utf-8").splitlines()

        title = slug
        for idx, line in enumerate(lines):
            if line.startswith("# "):
                title = line[2:].strip()
                del lines[idx]
                break

        body = normalise_headings(strip_local_paths(rewrite_links("\n".join(lines))))

        frontmatter = "\n".join([
            "---",
            f"title: {yaml_str(title)}",
            f"blurb: {yaml_str(blurb)}",
            f"order: {order}",
            "---",
            "",
        ])
        (out_dir / f"{slug}.md").write_text(frontmatter + body.lstrip("\n") + "\n", encoding="utf-8")

    print(f"md      : {len(PAGES)} pages -> src/content/replication/")


# ============================================================================
# 2. replication scripts
# ============================================================================


def sync_replication_scripts() -> list[dict]:
    out_dir = WEB / "src" / "replication-code"
    if out_dir.exists():
        for f in out_dir.glob("*.R"):
            f.unlink()
    out_dir.mkdir(parents=True, exist_ok=True)

    manifest = []
    seen = set()
    for f in sorted((ENH / "replication").rglob("*.R")):
        family = f.parent.name
        name = f.stem
        assert name not in seen, f"duplicate script name: {name}"
        seen.add(name)
        (out_dir / f.name).write_text(f.read_text(encoding="utf-8"), encoding="utf-8")
        manifest.append({"family": family, "name": name, "file": f.name})

    data_dir = WEB / "src" / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    (data_dir / "replication-scripts.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    print(f"scripts : {len(manifest)} -> src/replication-code/")
    return manifest


# ============================================================================
# 3. Rd parsing -- native reference pages, no pkgdown redirect
# ============================================================================

CMD_RE = re.compile(r"\\([A-Za-z]+)")

# Sections captured verbatim (one brace group) at the top level of an .Rd file.
# \section{title}{body} is repeatable and handled separately (two groups).
OUTER_SECTIONS = {
    "name", "alias", "title", "usage", "arguments", "value", "description",
    "details", "examples", "references", "seealso", "note", "format",
    "source", "author", "keyword", "docType",
}

# Rd markup that isn't real content -- roxygen2's header comment, e.g.
# "% Generated by roxygen2: do not edit by hand". R comments use '#', not
# '%', so a naive whole-file strip of '%'-led lines never touches real code.
_COMMENT_RE = re.compile(r"(?m)^%.*$\n?")


def strip_comments(text: str) -> str:
    return _COMMENT_RE.sub("", text)


def read_group(s: str, i: int) -> tuple[str, int]:
    """s[i] == '{'. Returns (content, index just past the matching '}')."""
    assert s[i] == "{"
    depth = 1
    j = i + 1
    n = len(s)
    while j < n:
        c = s[j]
        if c == "\\" and j + 1 < n and s[j + 1] in "{}":
            j += 2
            continue
        if c == "{":
            depth += 1
        elif c == "}":
            depth -= 1
            if depth == 0:
                return s[i + 1 : j], j + 1
        j += 1
    return s[i + 1 :], n  # unterminated -- be lenient, take the rest


def try_group(text: str, j: int) -> tuple[str | None, int]:
    if j < len(text) and text[j] == "{":
        return read_group(text, j)
    return None, j


# ---- inline parser: Rd markup -> a small JSON-serializable node tree -------
#
# node shapes:
#   {"t": "text", "v": str}
#   {"t": "code" | "strong" | "em", "c": [node, ...]}
#   {"t": "link", "topic": str|None, "external": bool, "href": str|None, "c": [node, ...]}
#   {"t": "eqn", "tex": str, "display": bool}
#   {"t": "list", "ordered": bool, "items": [{"name": [node,...]|None, "body": [node,...]}]}
#
# \eqn/\deqn keep the real LaTeX (their first group) rather than the ASCII
# fallback, so reference pages can render them with the same KaTeX component
# methodology.astro already uses.


def parse_inline(text: str) -> list[dict]:
    nodes: list[dict] = []
    buf: list[str] = []
    i, n = 0, len(text)

    def flush():
        if buf:
            nodes.append({"t": "text", "v": "".join(buf)})
            buf.clear()

    while i < n:
        c = text[i]
        if c != "\\":
            buf.append(c)
            i += 1
            continue
        m = CMD_RE.match(text, i)
        if not m:
            if i + 1 < n and text[i + 1] in "%_&#{}$^~\\":
                buf.append(text[i + 1])
                i += 2
            else:
                buf.append(c)
                i += 1
            continue
        name, j = m.group(1), m.end()
        node, j = parse_command(name, text, j)
        flush()
        if node is None:
            pass
        elif isinstance(node, list):
            nodes.extend(node)
        else:
            nodes.append(node)
        i = j
    flush()
    return nodes


def group_items(children: list[dict]) -> list[dict]:
    """Collects \\item nodes (braced or bare-marker) into {name, body} pairs."""
    items: list[dict] = []
    pending: dict | None = None
    for n in children:
        if n.get("t") == "_item_bare":
            if pending is not None:
                items.append(pending)
            pending = {"name": None, "body": []}
        elif n.get("t") == "item":
            if pending is not None:
                items.append(pending)
                pending = None
            items.append({"name": n["name"], "body": n["body"]})
        elif pending is not None:
            pending["body"].append(n)
        # else: stray content before the first \item (whitespace) -- drop
    if pending is not None:
        items.append(pending)
    return items


def parse_command(name: str, text: str, j: int):
    if name in ("eqn", "deqn"):
        g1, j = try_group(text, j)
        if g1 is None:
            return None, j
        _g2, j2 = try_group(text, j)  # ASCII fallback -- discarded, KaTeX renders g1
        return {"t": "eqn", "tex": g1.strip(), "display": name == "deqn"}, j2

    if name == "ifelse":
        _g1, j = try_group(text, j)
        _g2, j = try_group(text, j)
        g3, j = try_group(text, j)
        return (parse_inline(g3) if g3 is not None else []), j

    if name == "if":
        # \if{test}{content} -- HTML/latex-conditional content (a logo figure,
        # a lifecycle badge). We render as plain text, so -- like Rd2txt would
        # for a non-html/latex target -- the content is simply dropped.
        _g1, j = try_group(text, j)
        _g2, j = try_group(text, j)
        return None, j

    if name == "method":
        g1, j = try_group(text, j)
        _g2, j = try_group(text, j)
        return (parse_inline(g1) if g1 is not None else []), j

    if name == "item":
        g1, j = try_group(text, j)
        if g1 is None:
            return {"t": "_item_bare"}, j
        g2, j2 = try_group(text, j)
        if g2 is not None:
            return {"t": "item", "name": parse_inline(g1), "body": parse_inline(g2)}, j2
        return {"t": "item", "name": None, "body": parse_inline(g1)}, j

    if name == "link":
        # \link{topic} and \link[=dest]{text} stay internal. \link[pkg]{topic}
        # and \link[pkg:topic]{text} are cross-package (ggplot2::autoplot,
        # generics::augment, tibble::tibble, ...) -- no derivable href on this
        # site, so they render as plain code rather than a dead internal link.
        opt = None
        if j < len(text) and text[j] == "[":
            k = text.find("]", j)
            if k != -1:
                opt, j = text[j + 1 : k], k + 1
        g1, j = try_group(text, j)
        raw = (g1 or "").strip()
        disp = parse_inline(g1) if g1 is not None else []
        if opt is None:
            topic, external = (raw or None), False
        elif opt.startswith("="):
            topic, external = (opt[1:] or None), False
        else:
            topic, external = None, True
        return {"t": "link", "topic": topic, "external": external, "href": None,
                "c": disp or [{"t": "text", "v": raw}]}, j

    if name == "href":
        g1, j = try_group(text, j)
        g2, j = try_group(text, j)
        return {"t": "link", "topic": None, "external": True, "href": (g1 or "").strip(),
                "c": parse_inline(g2) if g2 is not None else []}, j

    if name in ("figure", "tabular"):
        _g1, j = try_group(text, j)
        _g2, j = try_group(text, j)
        return None, j  # not renderable statically -- drop (a lifecycle badge image)

    if name in ("code", "verb", "preformatted"):
        g1, j = try_group(text, j)
        return {"t": "code", "c": parse_inline(g1) if g1 is not None else []}, j

    if name in ("strong", "bold"):
        g1, j = try_group(text, j)
        return {"t": "strong", "c": parse_inline(g1) if g1 is not None else []}, j

    if name in ("emph", "dQuote", "sQuote"):
        g1, j = try_group(text, j)
        inner = parse_inline(g1) if g1 is not None else []
        if name == "dQuote":
            inner = [{"t": "text", "v": '"'}, *inner, {"t": "text", "v": '"'}]
        elif name == "sQuote":
            inner = [{"t": "text", "v": "'"}, *inner, {"t": "text", "v": "'"}]
        return {"t": "em", "c": inner}, j

    if name in ("url", "email", "doi"):
        g1, j = try_group(text, j)
        raw = (g1 or "").strip()
        href = {"email": f"mailto:{raw}", "doi": f"https://doi.org/{raw}"}.get(name, raw)
        return {"t": "link", "topic": None, "external": True, "href": href,
                "c": [{"t": "text", "v": raw}]}, j

    if name == "pkg":
        g1, j = try_group(text, j)
        return {"t": "code", "c": parse_inline(g1) if g1 is not None else []}, j

    if name in ("donttest", "dontrun", "dontshow"):
        g1, j = try_group(text, j)
        return (parse_inline(g1) if g1 is not None else []), j

    if name in ("itemize", "enumerate", "describe"):
        g1, j = try_group(text, j)
        children = parse_inline(g1) if g1 is not None else []
        return {"t": "list", "ordered": name == "enumerate", "items": group_items(children)}, j

    if name in ("keyword", "docType", "concept", "encoding"):
        _g1, j = try_group(text, j)
        return None, j

    if name == "cr":
        return {"t": "text", "v": "\n"}, j
    if name in ("dots", "ldots"):
        return {"t": "text", "v": "..."}, j
    if name == "R":
        return {"t": "text", "v": "R"}, j

    # Unknown command: transparently unwrap its group if it has one, else
    # emit the bare name so no stray backslash leaks into rendered prose.
    g1, j = try_group(text, j)
    if g1 is not None:
        return parse_inline(g1), j
    return {"t": "text", "v": name}, j


def flatten_text(nodes: list[dict]) -> str:
    out = []
    for n in nodes:
        t = n.get("t")
        if t == "text":
            out.append(n["v"])
        elif t in ("code", "strong", "em", "link"):
            out.append(flatten_text(n.get("c", [])))
        elif t == "eqn":
            out.append(n["tex"])
    return re.sub(r"\s+", " ", "".join(out)).strip()


def _normalise_paragraph(nodes: list[dict]) -> list[dict]:
    out = []
    for n in nodes:
        if n.get("t") == "text":
            v = re.sub(r"[ \t]*\n[ \t]*", " ", n["v"])
            v = re.sub(r" {2,}", " ", v)
            if v:
                out.append({"t": "text", "v": v})
        else:
            out.append(n)
    if out and out[0]["t"] == "text":
        out[0] = {"t": "text", "v": out[0]["v"].lstrip()}
    if out and out[-1]["t"] == "text":
        out[-1] = {"t": "text", "v": out[-1]["v"].rstrip()}
    return [n for n in out if not (n["t"] == "text" and n["v"] == "")]


def split_paragraphs(nodes: list[dict]) -> list[list[dict]]:
    """Splits a flat node list into paragraphs on blank-line boundaries."""
    paras: list[list[dict]] = []
    cur: list[dict] = []
    for n in nodes:
        if n.get("t") == "text" and re.search(r"\n[ \t]*\n", n["v"]):
            parts = re.split(r"\n[ \t]*\n+", n["v"])
            first, *rest = parts
            if first:
                cur.append({"t": "text", "v": first})
            paras.append(cur)
            cur = []
            for mid in rest[:-1]:
                paras.append([{"t": "text", "v": mid}])
            if rest and rest[-1]:
                cur = [{"t": "text", "v": rest[-1]}]
        else:
            cur.append(n)
    paras.append(cur)
    cleaned = [_normalise_paragraph(p) for p in paras]
    return [p for p in cleaned if p]


def parse_list_or_prose(nodes: list[dict]) -> dict:
    items = group_items(nodes)
    if items:
        return {
            "kind": "list",
            "items": [
                {"name": flatten_text(it["name"]) if it["name"] else "", "desc": it["body"]}
                for it in items
            ],
        }
    return {"kind": "prose", "paragraphs": split_paragraphs(nodes)}


def unwrap_donttest(raw: str) -> str:
    t = raw.strip("\n")
    m = re.match(r"^\s*\\(donttest|dontrun|dontshow)\{", t)
    if m:
        inner, end = read_group(t, m.end() - 1)
        if t[end:].strip() == "":
            return inner.strip("\n")
    return t.strip()


def clean_usage(raw: str) -> str:
    s = re.sub(r"\\method\{([^}]*)\}\{[^}]*\}", r"\1", raw)
    s = re.sub(r"\\l?dots\b", "...", s)
    return s.strip("\n")


def scan_top_level(text: str) -> list[tuple[str, str | None, str | None]]:
    text = strip_comments(text)
    out: list[tuple[str, str | None, str | None]] = []
    i, n = 0, len(text)
    while i < n:
        if text[i] != "\\":
            i += 1
            continue
        m = CMD_RE.match(text, i)
        if not m:
            i += 1
            continue
        name, j = m.group(1), m.end()
        if name == "section":
            g1, j = try_group(text, j)
            g2, j = try_group(text, j)
            out.append(("section", g1, g2))
            i = j
        elif name in OUTER_SECTIONS:
            g1, j = try_group(text, j)
            out.append((name, g1, None))
            i = j
        else:
            i = j  # some other backslash sequence between sections -- skip past it
    return out


def clean_identifier(raw: str) -> str:
    """\\name{}/\\alias{} bodies are raw, unparsed text -- Rd still escapes
    characters that would otherwise start a comment or a command
    (\\%>\\% for the pipe operator's own Rd page), so run them through the
    same escape handling as everywhere else rather than taking them literally."""
    return flatten_text(parse_inline(raw))


def quick_scan_name_alias(text: str) -> tuple[str, list[str]]:
    """Cheap first pass: just enough to build the topic/alias map."""
    top = scan_top_level(text)
    topic = next((clean_identifier(g) for n, g, _ in top if n == "name" and g), "")
    aliases = [clean_identifier(g) for n, g, _ in top if n == "alias" and g]
    return topic, aliases


def parse_rd_file(text: str) -> dict:
    top = scan_top_level(text)
    by_name: dict[str, list[str]] = {}
    sections: list[dict] = []
    for name, g1, g2 in top:
        if name == "section":
            sections.append({
                "title": flatten_text(parse_inline(g1)) if g1 else "",
                "body": split_paragraphs(parse_inline(g2)) if g2 else [],
            })
        elif g1 is not None:
            by_name.setdefault(name, []).append(g1)

    def prose(key: str):
        raw = by_name.get(key)
        return split_paragraphs(parse_inline(raw[0])) if raw else None

    def list_or_prose(key: str):
        raw = by_name.get(key)
        return parse_list_or_prose(parse_inline(raw[0])) if raw else None

    topic = clean_identifier(by_name.get("name", [""])[0])
    title_raw = by_name.get("title")
    usage_raw = by_name.get("usage")
    args_raw = by_name.get("arguments")
    examples_raw = by_name.get("examples")

    return {
        "topic": topic,
        "aliases": [clean_identifier(a) for a in by_name.get("alias", [])],
        "title": flatten_text(parse_inline(title_raw[0])) if title_raw else topic,
        "usage": clean_usage(usage_raw[0]) if usage_raw else None,
        "description": prose("description"),
        "details": prose("details"),
        "arguments": (
            [
                {"name": flatten_text(it["name"]) if it["name"] else "", "desc": it["body"]}
                for it in group_items(parse_inline(args_raw[0]))
            ]
            if args_raw
            else None
        ),
        "value": list_or_prose("value"),
        "format": list_or_prose("format"),
        "examples": unwrap_donttest(examples_raw[0]) if examples_raw else None,
        "references": prose("references"),
        "seealso": prose("seealso"),
        "source": prose("source"),
        "author": prose("author"),
        "sections": sections,
    }


# A canonical topic becomes a URL segment (/reference/<topic>). Almost every
# Rd \name is a plain identifier; the one exception in this corpus is
# pipe.Rd's own topic, the re-exported magrittr operator "%>%" -- not safe as
# a path segment, and not worth a dedicated page (keyword: internal, no
# content beyond "Pipe operator"). Anything that doesn't fit gets no page and
# no link, just plain code.
URL_SAFE_RE = re.compile(r"^[A-Za-z0-9_.-]+$")


def resolve_links_deep(obj, alias_topic: dict[str, str]) -> None:
    if isinstance(obj, dict):
        if obj.get("t") == "link" and obj.get("topic") and not obj.get("external"):
            resolved = alias_topic.get(obj["topic"])
            if resolved and URL_SAFE_RE.match(resolved):
                obj["topic"] = resolved
            else:
                obj["external"] = True
                obj["topic"] = None
        for v in obj.values():
            resolve_links_deep(v, alias_topic)
    elif isinstance(obj, list):
        for v in obj:
            resolve_links_deep(v, alias_topic)


def build_reference() -> None:
    rd_files = sorted((EXUBER / "man").glob("*.Rd"))
    rd_texts = {f: f.read_text(encoding="utf-8") for f in rd_files}

    alias_topic: dict[str, str] = {}
    topic_title: dict[str, str] = {}
    for f, text in rd_texts.items():
        topic, aliases = quick_scan_name_alias(text)
        for a in [topic, *aliases]:
            alias_topic[a] = topic

    topics: dict[str, dict] = {}
    for f, text in rd_texts.items():
        doc = parse_rd_file(text)
        topic_title[doc["topic"]] = doc["title"]
        if URL_SAFE_RE.match(doc["topic"]):
            topics[doc["topic"]] = doc
        # else: no dedicated page (see URL_SAFE_RE) -- still resolvable as a
        # title via topic_title, just not a link target.

    resolve_links_deep(topics, alias_topic)

    pk = yaml.safe_load((EXUBER / "_pkgdown.yml").read_text(encoding="utf-8"))

    groups = []
    for g in pk["reference"]:
        desc = re.sub(r"\s+", " ", (g.get("desc") or "")).strip()
        m = re.search(r"docs/enhancements/([a-z-]+)\.md", desc)
        replication = m.group(1) if m and m.group(1) in PUBLISHED else None
        desc = re.sub(r";?\s*see docs/enhancements/[a-z-]+\.md", "", desc).strip()

        contents = g.get("contents") or []
        entries = []
        for name in contents:
            topic = alias_topic.get(name, name)
            entries.append({"name": name, "topic": topic, "title": topic_title.get(topic, "")})
        groups.append({"title": g["title"], "desc": desc, "replication": replication, "entries": entries})

    namespace = (EXUBER / "NAMESPACE").read_text(encoding="utf-8").splitlines()
    exports = [
        re.sub(r'^export\(|\)$|"', "", line)
        for line in namespace
        if line.startswith("export(")
    ]
    grouped_names = {e["name"] for g in groups for e in g["entries"]}
    grouped_topics = {e["topic"] for g in groups for e in g["entries"]}
    grouped_aliases = {
        a
        for topic, doc in topics.items()
        if topic in grouped_topics
        for a in [topic, *doc["aliases"]]
    }
    ungrouped = sorted(set(exports) - grouped_names - grouped_aliases)

    data_dir = WEB / "src" / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    (data_dir / "reference.json").write_text(
        json.dumps({"r": groups, "ungrouped": ungrouped, "topics": topics}, indent=2),
        encoding="utf-8",
    )

    n_topics = sum(len(g["entries"]) for g in groups)
    print(f"reference: {len(groups)} groups, {n_topics} topics, {len(topics)} Rd files parsed -> src/data/reference.json")
    if ungrouped:
        print(f"  {len(ungrouped)} exports with no _pkgdown.yml reference entry:")
        print(f"    {', '.join(ungrouped)}")


def main() -> None:
    sync_markdown()
    sync_replication_scripts()
    build_reference()


if __name__ == "__main__":
    sys.exit(main())
