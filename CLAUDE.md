# website

Astro + Tailwind v4 public site (methodology, calculation-suite docs, the
replication record, a generated reference index), deployed at
[exuber.kvasilopoulos.com](https://exuber.kvasilopoulos.com/).

## Commands

```sh
npm run dev           # local dev server
npm run build         # production build
npm run preview       # serve the build locally
npm run check:links   # scripts/check-links.mjs
npm run sync          # scripts/sync_content.py -- see below
```

## Content sync: generated files, don't hand-edit

This repo doesn't contain its own source content for two sections — it
vendors it from the project root, which lives outside this git repo and
outside Netlify's build (Netlify only checks out `website/`):

- `docs/enhancements/*.md` (curated subset, the `PAGES` list at the top of
  `scripts/sync_content.py`) → `src/content/replication/*.md`, plus every
  `docs/enhancements/replication/**/*.R` → `src/replication-code/*.R`
- `exuber/_pkgdown.yml` + every `exuber/man/*.Rd` (hand-rolled Rd parser,
  not `tools::parse_Rd()`) → `src/data/reference.json` (usage, arguments,
  value, examples, seealso — full content, so `/reference/<topic>` renders
  natively instead of linking out to pkgdown)

**Edit the source** (`docs/enhancements/*.md` or `exuber/R/*.R`'s roxygen
comments), never `src/content/replication/*.md` or `src/data/reference.json`
directly — re-run `npm run sync` (needs PyYAML: `pip install pyyaml`) and
`npm run build && npm run check:links` after. Since the source lives
outside this repo, this is a manual step, not CI — nothing keeps the
vendored copy in sync automatically.

## Deploy mechanism

Netlify builds and deploys on every push to `main`, triggered by a GitHub
webhook — no separate deploy step, no staging gate. **Pushing to `main`
ships to production immediately.** Treat it accordingly: confirm before
pushing, and prefer a branch/PR if the change needs review first.
