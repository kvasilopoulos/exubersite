#!/usr/bin/env Rscript

# Vendors curated content from the project root into website/.
#
# website/ is its own git repo and Netlify builds it from GitHub, but
# docs/enhancements/ and exuber/ live at the project root, outside it --
# nothing there exists at Netlify build time. So this runs manually and its
# output is committed. No submodule, no CI step, no Netlify build change.
#
#   Rscript website/scripts/sync-content.R
#
# Three jobs:
#   1. curated docs/enhancements/*.md  -> src/content/replication/*.md
#   2. docs/enhancements/replication/**/*.R -> src/replication-code/*.R
#   3. exuber/_pkgdown.yml + man/*.Rd  -> src/data/reference.json

suppressWarnings(suppressMessages({
  stopifnot(requireNamespace("yaml", quietly = TRUE))
  stopifnot(requireNamespace("jsonlite", quietly = TRUE))
}))

# ---- paths ------------------------------------------------------------------

script_path <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) normalizePath(f) else file.path(getwd(), "website/scripts/sync-content.R")
})
web  <- dirname(dirname(script_path))
root <- dirname(web)
enh  <- file.path(root, "docs", "enhancements")

stopifnot(dir.exists(enh), dir.exists(file.path(root, "exuber", "man")))

# Curated set. Deliberately withheld: papers/ (52 copyrighted PDFs),
# status.md (conversational dump), README.md / SUMMARY.md (internal bundle
# log), open-research-directions.md and practitioner-guidance.md (backlog,
# not replication).
pages <- list(
  list(
    slug  = "volatility-robustness",
    blurb = "Tests robust to time-varying innovation variance: time-transformed, kernel-purged, WLS, sign-based, and stochastic-coefficient routes."
  ),
  list(
    slug  = "dating-and-root-inference",
    blurb = "Origination, collapse and recovery dates, plus confidence intervals on the explosive root itself."
  ),
  list(
    slug  = "monitoring",
    blurb = "Sequential and real-time detection: training-vs-monitoring orchestration, CUSUM families, and closed-form boundaries."
  ),
  list(
    slug  = "multivariate",
    blurb = "Panel and cross-series tests -- common bubbles, co-bubbles, and bubble contagion."
  ),
  list(
    slug  = "alternative-paradigms",
    blurb = "Non-ADF-family approaches, principally the quantile-based global test and its recursive monitoring extension."
  ),
  list(
    slug  = "simulation-dgps",
    blurb = "Data-generating processes for the axes the original sim_*() functions do not cover."
  ),
  list(
    slug  = "references",
    blurb = "Full bibliography behind the replication record, organised by methodological family."
  )
)
published <- vapply(pages, `[[`, "", "slug")

# ---- link rewriting ---------------------------------------------------------

# The source markdown links into papers/*.pdf, into sibling .md files, and into
# the replication/*.R scripts. Rewriting here rather than in a remark plugin
# keeps the published markdown inspectable in git and costs no dependency.
rewrite_one <- function(lnk) {
  text   <- sub("^\\[((?s).*)\\]\\([^)]*\\)$", "\\1", lnk, perl = TRUE)
  target <- sub("^\\[(?s).*\\]\\(([^)]*)\\)$", "\\1", lnk, perl = TRUE)

  if (grepl("^(https?:|mailto:|#)", target)) return(lnk)

  bare   <- sub("^(\\.\\./)+", "", target)
  anchor <- if (grepl("#", bare)) sub("^[^#]*", "", bare) else ""
  path   <- sub("#.*$", "", bare)

  if (path == "") return(lnk)
  if (grepl("^papers/", path)) return(text)                    # withheld
  if (grepl("\\.md$", path)) {
    slug <- sub("\\.md$", "", basename(path))
    if (slug %in% published) return(sprintf("[%s](/replication/%s%s)", text, slug, anchor))
    return(text)                                               # unpublished
  }
  if (grepl("^replication/.*\\.R$", path)) {
    return(sprintf("[%s](#script-%s)", text, sub("\\.R$", "", basename(path))))
  }
  text                                                         # any other local path
}

rewrite_links <- function(txt) {
  m <- gregexpr("\\[[^]]*\\]\\([^)[:space:]]+\\)", txt, perl = TRUE)
  found <- regmatches(txt, m)[[1]]
  if (length(found)) {
    regmatches(txt, m) <- list(vapply(found, rewrite_one, "", USE.NAMES = FALSE))
  }
  txt
}

# Bare prose references to withheld or project-root paths -- "Local copy:
# papers/x/y.pdf.", "open — papers/x/y.pdf", "`docs/enhancements/replication/
# fam/script.R`". Internal bookkeeping about where a local file sits; it means
# nothing to a visitor and points at files the site does not ship.
strip_local_paths <- function(txt) {
  txt <- gsub(
    paste0("(?:[,;]?\\s*(?:[Ll]ocal\\s+cop(?:y|ies):?|[Ll]ocal\\s+access:?|",
           "[Cc]opy:?|[Aa]ccess:?)?\\s*(?:[-–—]\\s*)?)",
           "papers/\\S+?\\.(?:pdf|txt)[,.]?"),
    "", txt, perl = TRUE
  )
  txt <- gsub("docs/enhancements/replication/[a-z-]+/([A-Za-z0-9_]+)\\.R",
              "\\1.R", txt, perl = TRUE)
  txt <- gsub("docs/enhancements/([a-z-]+)\\.md", "/replication/\\1", txt, perl = TRUE)
  txt <- gsub("\\s*PDFs: `papers/`\\.", "", txt, perl = TRUE)
  txt <- gsub("`papers/`", "the paper library", txt, fixed = TRUE)
  gsub("[ \t]+\n", "\n", txt)
}

# references.md titles its sections "## Multivariate — [multivariate.md](...)".
# The trailing link pushes the heading slug to `multivariate--multivariatemd`,
# which is not what the sibling files link to. Move it onto its own line so the
# slug is the section name; run scripts/check-links.mjs after a build to catch
# any anchor that still does not resolve.
normalise_headings <- function(txt) {
  txt <- gsub(
    paste0("(?m)^(#{2,3} [^\n\\[]+?)\\s+[—–-]\\s+",
           "(?:\\[[^\\]]+\\]\\([^)]*\\)|[a-z0-9-]+\\.md)\\s*$"),
    "\\1", txt, perl = TRUE
  )
  # Stale in the source: alternative-paradigms.md still points at the anchor
  # references.md used before it was reorganised out of tier ordering.
  gsub("#tier-7--quantile-based-detection", "#alternative-paradigms", txt, fixed = TRUE)
}

yaml_str <- function(x) paste0('"', gsub('"', '\\\\"', x), '"')

# ---- 1. markdown ------------------------------------------------------------

out_md <- file.path(web, "src", "content", "replication")
unlink(out_md, recursive = TRUE)
dir.create(out_md, recursive = TRUE, showWarnings = FALSE)

for (i in seq_along(pages)) {
  p   <- pages[[i]]
  src <- file.path(enh, paste0(p$slug, ".md"))
  stopifnot(file.exists(src))
  lines <- readLines(src, warn = FALSE, encoding = "UTF-8")

  # Take the file's own H1 as the page title and drop it -- the layout renders
  # the title, so leaving it in would print it twice.
  h1 <- which(grepl("^# ", lines))[1]
  title <- if (!is.na(h1)) sub("^# ", "", lines[h1]) else p$slug
  if (!is.na(h1)) lines <- lines[-h1]

  body <- normalise_headings(strip_local_paths(rewrite_links(paste(lines, collapse = "\n"))))

  writeLines(
    c(
      "---",
      paste0("title: ", yaml_str(trimws(title))),
      paste0("blurb: ", yaml_str(p$blurb)),
      paste0("order: ", i),
      "---",
      "",
      trimws(body, which = "left")
    ),
    file.path(out_md, paste0(p$slug, ".md")),
    useBytes = TRUE
  )
}
cat(sprintf("md      : %d pages -> src/content/replication/\n", length(pages)))

# ---- 2. replication scripts -------------------------------------------------

out_r <- file.path(web, "src", "replication-code")
unlink(out_r, recursive = TRUE)
dir.create(out_r, recursive = TRUE, showWarnings = FALSE)

r_src <- list.files(file.path(enh, "replication"), pattern = "\\.R$",
                    recursive = TRUE, full.names = TRUE)
manifest <- lapply(r_src, function(f) {
  family <- basename(dirname(f))
  name   <- sub("\\.R$", "", basename(f))
  file.copy(f, file.path(out_r, basename(f)), overwrite = TRUE)
  list(family = family, name = name, file = basename(f))
})
stopifnot(!anyDuplicated(vapply(manifest, `[[`, "", "name")))

dir.create(file.path(web, "src", "data"), recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(manifest, file.path(web, "src", "data", "replication-scripts.json"),
                     auto_unbox = TRUE, pretty = TRUE)
cat(sprintf("scripts : %d -> src/replication-code/\n", length(manifest)))

# ---- 3. reference index -----------------------------------------------------

pk <- yaml::read_yaml(file.path(root, "exuber", "_pkgdown.yml"))

rd_files <- list.files(file.path(root, "exuber", "man"), pattern = "\\.Rd$", full.names = TRUE)
rd <- lapply(rd_files, function(f) {
  doc  <- tools::parse_Rd(f, encoding = "UTF-8")
  tags <- vapply(doc, function(x) attr(x, "Rd_tag"), "")
  flat <- function(i) gsub("\\s+", " ", trimws(paste(unlist(doc[[i]]), collapse = "")))
  ti   <- which(tags == "\\title")
  list(
    topic   = sub("\\.Rd$", "", basename(f)),
    title   = if (length(ti)) flat(ti[1]) else "",
    aliases = vapply(which(tags == "\\alias"), flat, "")
  )
})

# alias -> topic, so `index` resolves to index-rd.Rd and S3 methods resolve too
alias_topic <- unlist(lapply(rd, function(x) setNames(rep(x$topic, length(x$aliases)), x$aliases)))
topic_title <- setNames(vapply(rd, `[[`, "", "title"), vapply(rd, `[[`, "", "topic"))

groups <- lapply(pk$reference, function(g) {
  desc  <- gsub("\\s+", " ", trimws(g$desc %||% ""))
  # Several groups already point at docs/enhancements/<slug>.md in their desc;
  # that mapping is 1:1 with the replication pages, so surface it as a link.
  slug  <- sub(".*docs/enhancements/([a-z-]+)\\.md.*", "\\1", desc)
  slug  <- if (slug != desc && slug %in% published) slug else NA_character_
  desc  <- trimws(sub("; see docs/enhancements/[a-z-]+\\.md", "", desc))

  entries <- lapply(unlist(g$contents), function(nm) {
    topic <- if (!is.na(alias_topic[nm])) unname(alias_topic[nm]) else nm
    ttl   <- unname(topic_title[topic])
    list(name = nm, topic = topic, title = if (is.na(ttl)) "" else ttl)
  })
  list(title = g$title, desc = desc, replication = slug, entries = entries)
})

exports <- local({
  ns <- readLines(file.path(root, "exuber", "NAMESPACE"), warn = FALSE)
  e  <- grep("^export\\(", ns, value = TRUE)
  gsub('^export\\(|\\)$|"', "", e)
})
grouped <- unlist(lapply(pk$reference, function(g) unlist(g$contents)))
grouped_aliases <- unlist(lapply(rd[vapply(rd, function(x) x$topic %in%
  unname(ifelse(is.na(alias_topic[grouped]), grouped, alias_topic[grouped])), TRUE)],
  `[[`, "aliases"))
ungrouped <- sort(setdiff(exports, union(grouped, grouped_aliases)))

jsonlite::write_json(
  list(r = groups, ungrouped = ungrouped),
  file.path(web, "src", "data", "reference.json"),
  auto_unbox = TRUE, pretty = TRUE
)
cat(sprintf("reference: %d groups, %d topics -> src/data/reference.json\n",
            length(groups), length(grouped)))
if (length(ungrouped)) {
  cat(sprintf("  %d exports with no _pkgdown.yml reference entry:\n    %s\n",
              length(ungrouped), paste(ungrouped, collapse = ", ")))
}
