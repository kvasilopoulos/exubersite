# website

An Astro and Tailwind v4 site with the methodology, the calculation-suite
docs, the replication record and a generated reference index. It is deployed
at [exuber.kvasilopoulos.com](https://exuber.kvasilopoulos.com/).

## Commands

```sh
npm run dev           # local dev server
npm run build         # production build
npm run preview       # serve the build locally
npm run check:links   # scripts/check-links.mjs
npm run sync          # scripts/sync_content.py, see below
```

## Content sync: generated files, do not hand-edit

Two sections of the site have no source content in this repo. Their content is
vendored from the project root, which sits outside this git repo and outside
the Netlify build (Netlify only checks out `website/`):

- The curated subset of `docs/*.md` (the `PAGES` list at the top of
  `scripts/sync_content.py`) goes to `src/content/replication/*.md`. Every
  `docs/replication/**/*.R` goes to `src/replication-code/*.R`.
- The "What's actually implemented" table in `docs/README.md` goes to
  `src/data/replication-table.json`, the flat per-method index on
  `/replication`. Its `Family` column must be a published page slug.
- `exuber/_pkgdown.yml` and every `exuber/man/*.Rd` go to
  `src/data/reference.json` through a hand-written Rd parser (not
  `tools::parse_Rd()`). The file holds the full usage, arguments, value,
  examples and seealso, so `/reference/<topic>` renders on this site and does
  not link out to pkgdown.

Edit the source, which is `docs/*.md` or the roxygen comments in
`exuber/R/*.R`. Never edit `src/content/replication/*.md` or
`src/data/reference.json` directly. After a source change, re-run
`npm run sync` (it needs PyYAML: `pip install pyyaml`), then
`npm run build && npm run check:links`. The source lives outside this repo, so
this is a manual step and not part of CI. Nothing updates the vendored copy
automatically.

## Deploy mechanism

Netlify builds and deploys on every push to `main`, triggered by a GitHub
webhook. There is no separate deploy step and no staging gate, so pushing to
`main` ships to production immediately. Confirm before pushing, and prefer a
branch and pull request if the change needs review first.

## Writing style (all user-facing text)

Applies to READMEs, vignettes, the website, `docs/`, NEWS/CHANGELOG,
roxygen and docstrings, and any prose a reader sees. Code comments and
CLAUDE.md files follow it too.

**Voice.** An applied economist writing for colleagues who also want
ordinary readers to be able to run the test. Precise, sober, a little
plain-spoken. Define a term at first use (what "explosive" means, what a
critical value is for) and give the idea in words before the formula.

**Rewrite, do not substitute.** Swapping an em dash for a comma, colon or
hyphen keeps the machine-written rhythm and is not acceptable. If a
sentence needed a dash, it was carrying two thoughts: split it into two
sentences, or fold the aside into the grammar (a relative clause, a
parenthesis only for a true aside, or a separate sentence). No U+2014 and
no spaced hyphen standing in for one. En dashes stay for numeric ranges
and joint names (Phillips–Shi–Yu).

**Patterns to remove at the sentence level.**
- Fragments stacked for effect, and "X, not Y" or "not just X, but Y"
  framings. State the claim directly.
- Triplets used for rhythm, and sentences that announce what they are about
  to say ("Importantly,", "It is worth noting that", "In essence").
- Telegraphic notes (dropped articles, arrows, semicolon chains, "confirmed,
  zero new code"). Write full sentences with a subject and a verb.
- Status-report voice: "genuinely", "confirmed", "now done", "picked
  clean", "the most topical candidate". Say what is true and give the date
  if it matters.
- Marketing and filler words: seamlessly, robust (unless a statistical
  sense is stated), leverage, delve, comprehensive, powerful, crucial,
  landscape, journey, "under the hood", "a rich set of".
- Hedge stacks, and bold used as emphasis inside running prose.
- Self-reference to the writing process ("this resolves the question this
  file flagged", "an earlier pass"). Keep history in dated notes, not in
  the body of explanations.

**Do keep.** Formulas, numbers, citations, function names and every fact.
This is a change of language, not of content. Vary sentence length. Prefer
"we" or the imperative to the passive, and say what a function does and
when to use it before how it works.
