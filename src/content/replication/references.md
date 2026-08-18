---
title: "References — full bibliography"
blurb: "Full bibliography behind the replication record, organised by methodological family."
order: 7
---
Every paper touched by this project, in one place, organized to match the
taxonomy files (`volatility-robustness.md`, `dating-and-root-inference.md`,
etc. — see `README.md`). Per-item detail, formulas, and exact-number
verification live in those files; this is the index. "Access" reflects
what was actually found during research, not a general claim about the
paper's OA status.

the paper library holds local PDF copies. As of 2026-08-09, **all 47 papers have a
local copy** (linked inline). 31 came via open-access routes, 14 via
institutional access once the session was routed through a UK academic
network (Jisc/Lancaster) — see "What changed between the first and second
pass" below for the mechanics of each publisher. The last 2 (HLS 2017,
HLW 2020, both Elsevier) resisted every automated route, including
institutional-network Playwright — see "The two that automation couldn't
get" below. 4 more were added 2026-08-11 specifically for
[simulation-dgps.md](/replication/simulation-dgps)'s DGP survey (bringing the total
to 52), all open access.

## Foundational (pre-existing in exuber)

| Citation | Access | Role |
|---|---|---|
| Phillips, P.C.B., Wu, Y. & Yu, J. (2011). "Explosive Behavior in the 1990s Nasdaq: When Did Exuberance Escalate Asset Values?" *International Economic Review*, 52(1), 201-226. Working paper: Cowles Foundation DP 1699. | open | PWY sup-ADF — `radf()`'s `sadf`/`badf` |
| Phillips, P.C.B., Shi, S. & Yu, J. (2015a). "Testing for Multiple Bubbles: Historical Episodes of Exuberance and Collapse in the S&P 500." *International Economic Review*, 56(4), 1043-1078. | open | PSY/GSADF — `radf()`'s `gsadf` |
| Phillips, P.C.B., Shi, S. & Yu, J. (2015b). "Testing for Multiple Bubbles: Limit Theory of Real-Time Detectors." *International Economic Review*, 56(4), 1079-1134. Working paper: Cowles Foundation DP 1915. | open | BSADF sequence + `datestamp()` |
| Pavlidis et al. (2016) (cited by exuber's existing code, not re-verified) | not chased | sieve bootstrap already in `radf_sb.R`/`radf_sb_cv()` |

## Volatility robustness
| # | Citation | Access | Status |
|---|---|---|---|
| 1 | Harvey, D.I., Leybourne, S.J. & Zu, Y. (2020). "Sign-based unit root tests for explosive financial bubbles in the presence of deterministically time-varying volatility." *Econometric Theory*, 36(1), 122-169. `doi:10.1017/S0266466619000057` | open (institutional access) via Cambridge Core. Its 2025 OBES level-shift extension also recovered | todo |
| 2 | Kurozumi, E., Skrobotov, A. & Tsarev, A. (2024). "Time-Transformed Test for Bubbles under Non-stationary Volatility." *J. Financial Econometrics*. `doi:10.1093/jjfinec/nbae026` | open (arXiv:2012.13937) | **done** |
| 3 | Harvey, D.I., Leybourne, S.J. & Zu, Y. (2019). "Testing explosive bubbles with time-varying volatility." *Econometric Reviews*, 38(10), 1131-1151. | open (Nottingham Granger Centre DP 18/05) | **done** |
| 4 | Harvey, D.I., Leybourne, S.J., Taylor, A.M.R. & Zu, Y. (2024/2025). "A new heteroskedasticity-robust test for explosive bubbles." *JTSA*, 46(5), 846-866. `doi:10.1111/jtsa.12784` | open (CC-BY) | **done** (with-intercept variant) |
| 5 | Hafner, C.M. (2020). "Testing for Bubbles in Cryptocurrencies with Time-Varying Volatility." *J. Financial Econometrics*, 18(2), 233-249. | open (EconStor/IRTG 1792 DP 2018-005, SSRN itself blocks automated fetch) | evaluated, not implemented |
| 6 | Pedersen, T.Q. & Montes Schütte, E.C. (2020). "Testing for Explosive Bubbles in the Presence of Autocorrelated Innovations." *J. Empirical Finance*, 58, 207-225. | open (Aarhus CREATES RP 2017-9) | evaluated — mostly already covered by `radf_sb_cv()` |
| 7 | Kurozumi, E. & Nishi, M. (2025). "Testing for a bubble with a stochastically varying explosive coefficient." *JTSA*, 46(5), 945-965. | **open (Open Access article)**; it wasn't paywalled at all, Wiley just requires JS to serve it, which plain `curl` can't do | todo, newly surfaced |
| 8 | Sarkar, A. & Wells, M.T. (2026). "Is There an AI Bubble? Robust Date-Stamping for Periods of Exuberance." arXiv:2604.12062 (preprint). Theory: Sarkar & Wells, "Double Local-to-Unity," arXiv:2512.06823. | open + | todo, preprint |
| 9 | Monschang, V. & Wilfling, B. (2021). "Sup-ADF-style bubble-detection methods under test." *Empirical Economics*, 61, 145-172. `doi:10.1007/s00181-020-01859-7` | open (CQE WP 78/2019) | new (2026-08-11) — TGARCH leverage + lognormal-mixture bubble DGP, see [simulation-dgps.md](/replication/simulation-dgps#12-tgarch11-with-leverage-effect) |
| 10 | Richter, S., Wang, W. & Wu, W.B. (2023, orig. 2018). "A supreme test for periodic explosive GARCH." *Econometrics* (MDPI); arXiv:1812.03475. | open | new (2026-08-11) — volatility-bubble (not price-bubble) GARCH, out of scope for `radf()`, see [simulation-dgps.md](/replication/simulation-dgps#adjacent-not-a-price-level-dgp) |

## Dating and root inference
| Citation | Access | Notes |
|---|---|---|
| Harvey, D.I., Leybourne, S.J. & Sollis, R. (2017). "Improving the accuracy of asset price bubble start and end date estimators." *J. Empirical Finance*, 40, 121-138. `doi:10.1016/j.jempfin.2016.11.001` | open (supplied by user; automated fetch never got past Elsevier/Nottingham's Cloudflare gate — see roundup below) | HLS |
| Harvey, D.I., Leybourne, S.J. & Whitehouse, E.J. (2020). "Date-stamping multiple bubble regimes." *J. Empirical Finance*, 58, 226-246. `doi:10.1016/j.jempfin.2020.06.004` | open (supplied by user; same blocker as HLS) | HLW |
| Pang, T., Du, L. & Chong, T.T.L. (2021). "Estimating multiple breaks in nonstationary autoregressive models." *J. Econometrics*, 221(1), 277-311. | open (MPRA 92074) | PDC |
| Kurozumi, E. & Skrobotov, A. (2023). "On the asymptotic behavior of bubble date estimators." *JTSA*, 44(4), 359-373. | open (arXiv:2110.04500) | KS |
| Kurozumi, E. & Skrobotov, A. (2023). "Improving the accuracy of bubble date estimators under time-varying volatility." arXiv:2306.02977. | open | newly found — two-step WLS-based dating estimator |
| Kejriwal, M., Nguyen, L. & Perron, P. (2025). "An Improved Procedure for Retrospectively Dating the Emergence and Collapse of Bubbles." *JTSA*, 46(5). `doi:10.1111/jtsa.12810` | open (institutional access) | newly surfaced |
| Kurozumi, E. & Skrobotov, A. (2025). "Confidence Sets for the Emergence, Collapse, and Recovery Dates of a Bubble." arXiv:2511.16172 | open | newly surfaced |
| Phillips, P. C. B., & Magdalinos, T. (2007). "Limit theory for moderate deviations from a unit root." *J. Econometrics*, 136(1), 115-130. Working paper: Cowles Foundation DP 1471 (July 2004). | **open** — Theorem 4.3 verified by PNG page-render 2026-08-09 | source-verified, not yet coded |
| Guo, G., Sun, Y. & Wang, S. (2019). "Testing for moderate explosiveness." *The Econometrics Journal*, 22(3), 279-303. | open (institutional access); `rootstamp()` was implemented from Skrobotov's secondary restatement before this was found — worth cross-checking the implementation against the primary source now that it's available | **done** — `rootstamp()` |
| Phillips, P.C.B., Magdalinos, T. & Giraitis, L. (2010). "Smoothing local-to-moderate unit root theory." *J. Econometrics*, 158(2), 274-279. Working paper: Cowles Foundation DP 1659. | open | background theory only |
| Phillips, P.C.B. & Shi, S. (2014). "Financial Bubble Implosion." Working paper: Cowles Foundation DP 1967. Published as "Financial Bubble Implosion and Reverse Regression," *Econometric Theory*. | open | newly surfaced — reverse-BSADF recovery dating |

## Monitoring
| Citation | Access | Notes |
|---|---|---|
| Homm, U. & Breitung, J. (2012). "Testing for speculative bubbles in stock markets: a comparison of alternative methods." *J. Financial Econometrics*, 10(1), 198-231. `doi:10.1093/jjfinec/nbr009` | open (institutional access) | formulas now directly readable, not yet re-verified against this copy |
| Astill, S., Harvey, D.I., Leybourne, S.J., Taylor, A.M.R. & Zu, Y. (2021/2023). "CUSUM-Based Monitoring for Explosive Episodes in Financial Data in the Presence of Time-Varying Volatility." *J. Financial Econometrics*, 21(1), 187-227. `doi:10.1093/jjfinec/nbab009` | open (institutional access) | AHLTZ; primary source now available, not yet re-verified against this copy |
| Astill, S., Harvey, D.I., Leybourne, S.J., Sollis, R. & Taylor, A.M.R. (2018). "Real-Time Monitoring for Explosive Financial Bubbles." *JTSA*, 39, 863-891. | open (institutional access) | AHLST |
| Whitehouse, E.J., Harvey, D.I. & Leybourne, S.J. (2025). "Real-time monitoring procedures for early detection of bubbles." *Intl J. Forecasting*, 41(3), 1260-1277. `doi:10.1016/j.ijforecast.2024.12.005` | **open (CC-BY)** | supplied verified AHLST decision rule + FPR formula |
| Kurozumi, E. (2020). "Asymptotic properties of bubble monitoring tests." *Econometric Reviews*, 39(5), 510-538. `doi:10.1080/07474938.2019.1697086` | open (institutional access) | now primary-source readable, was abstract-level only before |
| Kurozumi, E. (2021). "Asymptotic Behavior of Delay Times of Bubble Monitoring Tests." *JTSA*, 42(3), 314-337. `doi:10.1111/jtsa.12569` | open (institutional access) | now primary-source readable, was abstract-level only before |
| Horváth, L. & Trapani, L. (2026). "Real-time monitoring with RCA models." *Econometric Theory*, 42, 514-547. Working paper arXiv:2312.11710. | open | |
| Astill, S., Taylor, A.M.R. & Zu, Y. (2026, forthcoming). "Covariate Augmented CUSUM Bubble Monitoring Procedures." *Econometric Theory*. Working paper: Essex Finance Centre WP 94. | open | most useful single doc found for this file |
| Breitung, J. & Diegel, M. (2025). "Sequential Detector Statistics for Speculative Bubbles." *JTSA*, 46(5). | open (via EconStor) | newly surfaced — direct match |

## Multivariate
| Citation | Access | Notes |
|---|---|---|
| Chen, Y., Phillips, P.C.B. & Shi, S. (2020/2023). "Common Bubble Detection in Large Dimensional Financial Systems." *J. Financial Econometrics*, 21(4), 989-1063. Working paper: Cowles Foundation DP 2251. | **open** | **done** — `radf_common()` |
| Evripidou, A.C., Harvey, D.I., Leybourne, S.J. & Sollis, R. (2022). "Testing for Co-explosive Behaviour in Financial Time Series." *OBES*, 84(3), 624-650. `doi:10.1111/obes.12487` | open (institutional access) | evaluated, blocked — now unblocked, worth re-reading from primary source |
| Phillips, P.C.B. & Yu, J. (2011). "Dating the Timeline of Financial Bubbles During the Subprime Crisis." *Quantitative Economics*, 2(3), 455-491. `doi:10.3982/QE82` | open (Cowles DP 1770) | evaluated — no new code needed |
| Greenaway-McGrevy, R. & Phillips, P.C.B. (2015/2016). "Hot Property in New Zealand: Empirical Evidence of Housing Bubbles in the Metropolitan Centres." *NZ Economic Papers*, 50(1), 88-113. Working paper: Cowles Foundation DP 2004. | open | evaluated, not implemented |

## Open research directions
| Citation | Access | Notes |
|---|---|---|
| Skrobotov, A. (2023). "Testing for explosive bubbles: a review." *Dependence Modeling*, 11(1), 1-26. | open (arXiv:2207.08249) | source of this whole project |
| Lui, Y.L., Phillips, P.C.B. & Yu, J. (2024). "Robust Testing for Explosive Behavior with Strongly Dependent Errors." *J. Econometrics*, 238(2), 105626. Working paper: Cowles Foundation DP 2350. | **open** | evaluated, not implemented |
| Qian, J. & Su, L. (2016). "Shrinkage Estimation of Regression Models With Multiple Structural Changes." *Econometric Theory*, 32(6), 1376-1433. `doi:10.1017/S0266466615000237` | open (institutional access) via Cambridge Core; not itself a bubble-specific application | no bubble-specific LASSO precedent found |

## Alternative paradigms
| Citation | Access | Notes |
|---|---|---|
| Pavlidis, E.G. (2025). "Bubbles and crashes: A tale of quantiles." *JTSA*, 46(5), 884-907. | open (Lancaster ePrints) | quantile-based |
| Wu, R., Shi, S. & Wu, J. (2025). "Quantile analysis for financial bubble detection and surveillance." *JTSA*, 46(5), 908-931. | open (institutional access); SSRN itself is still blocked, this came via Wiley instead | quantile-based |
| Blasques, F., Koopman, S.J., Mingoli, G. & Telg, S. (2025). "A Novel Test for the Presence of Local Explosive Dynamics." *JTSA*, 46(5), 966-980. `doi:10.1111/jtsa.70001` | open (Tinbergen DP 24-036/III) | noncausal |
| Gourieroux, C. & Jasiak, J. (2025). "A Stochastic Tree for Bubble Asset Modelling and Pricing." *JTSA*, 46(5), 932-944. | open (institutional access, and turns out to be Open Access itself) | out of scope — pricing, not testing |
| Bhandari, A. "Rational Bubbles at the Spectral Edge." arXiv:2607.03933. | open | out of scope — different paradigm entirely |
| Chan, J.C.C. & Santi, C. (2021). "Speculative Bubbles in Present-Value Models: A Bayesian Markov-Switching State Space Approach." *J. Economic Dynamics and Control*, 127, 104101. | open (author's site) | new (2026-08-11) — Markov-switching bubble/collapse timing, see [simulation-dgps.md](/replication/simulation-dgps#14-markov-switching-present-value-bubble) |
| Chen, H., Chen, L., Huang, D., Li, Y. & Zhang, Z. (2026). "Technology Fundamentals and False Bubble Detection: Evidence from Dot-Com and AI Episodes." arXiv:2604.25826. | open | new (2026-08-11) — deterministic false-bubble null, see [simulation-dgps.md](/replication/simulation-dgps#15-deterministic-technology-adoption-false-bubble-null) |

## Practitioner guidance
| Citation | Access |
|---|---|
| Phillips, P.C.B., Shi, S. & Yu, J. (2014). "Specification Sensitivity in Right-Tailed Unit Root Testing for Explosive Behaviour." *OBES*, 76(3), 315-333. | open |
| Phillips, P.C.B. & Shi, S. (2019). "Detecting Financial Collapse and Ballooning Sovereign Risk." *OBES*, 81(6), 1336-1361. Working paper: Cowles Foundation DP 2110. | open (institutional access); SSRN and the SMU repository copy are both still blocked, this came via Wiley instead |
| Basele, R.B., Phillips, P.C.B. & Shi, S. (2025). "Speculative Bubbles in the Recent AI Boom: Nasdaq and the Magnificent Seven." *JTSA*, 46(5). | open |

## Context / meta

- Harvey, D.I. & Leybourne, S.J. (2025). Guest editors' introduction, *JTSA* special issue "Recent Developments in Time-Series Methods for Detecting Bubbles and Crashes," 46(5). `doi:10.1111/jtsa.70003` — paywalled itself, but its existence is why so much of this project got a 2025-26 refresh; full issue ToC recovered via RePEc/IDEAS, not the paywalled editorial text.
- Hu, Y. et al. (2023). "A review of Phillips-type right-tailed unit root bubble detection tests." *J. Economic Surveys*, 37(1), 141-158. — a second review, alongside Skrobotov (2023), not yet used as a source for this project but flagged as a candidate second pass.

## The two that automation couldn't get: HLS and HLW

**Harvey, Leybourne & Sollis (2017)** and **Harvey, Leybourne & Whitehouse
(2020)** (both *J. Empirical Finance*, in
[dating-and-root-inference.md](/replication/dating-and-root-inference)) were the
last two holdouts, and both are now in the paper library anyway — supplied
directly by the user, not fetched by any tool here. Worth recording why
the automated routes failed, since it's the one clean counterexample in an
otherwise near-total recovery: both are Elsevier (ScienceDirect), and both
have an open Green-OA copy per Unpaywall at
`nottingham-repository.worktribe.com` — but that repository 403s every
request behind a Cloudflare "Just a moment…" challenge, and it held even
against a real Playwright browser session routed through the same academic
network that got everything else through. ScienceDirect's own article
pages hit the identical wall. Unlike the Nottingham *Granger Centre*
old-style hosting (which just serves files directly, no gate at all), the
newer institutional-repository platform apparently runs a stricter
bot-check than Wiley, OUP, Taylor & Francis, or Cambridge Core did — none
of which blocked the same browser session. The Wayback Machine has only
ever crawled the HTML landing pages for these two, never the underlying
PDF asset (checked via the CDX API), so that fallback didn't apply either.

## What changed between the first and second pass

The first pass (plain `curl`, no institutional network) found 31 of 47
papers and left 16 "not found." The second pass, run once the session was
routed through Jisc/Lancaster's academic network, recovered 14 of those 16
by the same DOI, via:

- **Cambridge Core**: worked immediately with plain `curl` once
  IP-recognized — no JS, no browser needed (2/2 recovered this way).
- **Wiley** (`onlinelibrary.wiley.com`): the site is a JS-rendered SPA
  regardless of subscription status — `curl` gets a 60KB Angular shell no
  matter what. Needed a real browser (Playwright) to load the page, then
  `page.evaluate(() => fetch(...))` from inside that page's JS context to
  pull the actual PDF bytes as base64 (Chrome's built-in PDF viewer
  intercepts a plain navigation before the raw bytes are capturable any
  other way) (8/8 recovered this way — one of the eight, Kurozumi & Nishi
  2025, turned out to be genuinely Open Access, not paywalled at all).
- **Oxford University Press** (`academic.oup.com`) and **Taylor & Francis**
  (`tandfonline.com`): also needed a real browser — first attempts hit the
  same Cloudflare challenge as Elsevier/Nottingham, but subsequent
  attempts against the *same* URLs got through cleanly, suggesting a
  probabilistic/rate-limited challenge rather than a hard block. PDFs were
  then on a different subdomain (`watermark02.silverchair.com`, OUP's PDF
  host) — same `fetch()`-from-page-context trick worked there too (4/4
  recovered this way, after retries).
- **Elsevier** (`sciencedirect.com`) and the **Nottingham repository**:
  never got past the Cloudflare challenge, browser or not — see above.

Net effect: this is less "here are more papers" and more "many formulas
that were previously sourced secondhand (via restatements, abstracts, or
Skrobotov's review) can now be verified against the actual primary text"
— several rows above are flagged "worth re-reading from primary source now
that it's available," which is real follow-up work, not just a
citation-list improvement.

## Recurring access routes (useful for future passes)

- **On an institutional network, prefer the publisher's own page over
  hunting for a working-paper mirror.** Once IP-recognized, Wiley/OUP/T&F/
  Cambridge Core all served full text directly — no need for the
  Cowles/arXiv/MPRA/EconStor detective work below, which matters mainly
  when *not* on such a network, or for the handful of publishers (Elsevier,
  and repositories fronted by the newer Cloudflare-grade bot-checks) that
  block automated access independent of subscription status.
- **Wiley** (`onlinelibrary.wiley.com`): needs a real browser — it's a
  JS-rendered SPA for everyone, subscriber or not. `doi/pdfdirect/{doi}`
  is the endpoint; fetch it via `page.evaluate(() => fetch(...))` from
  within an already-loaded Wiley page (same-origin, so no CORS issue),
  base64-encode the bytes in-page, and decode outside the browser —
  navigating to the URL directly instead just hands the bytes to Chrome's
  built-in PDF viewer extension, which intercepts them before they're
  capturable as a plain network response body.
- **OUP** (`academic.oup.com`) and its PDF host **Silverchair**
  (`watermark02.silverchair.com`): same `fetch()`-from-page-context
  technique works once past OUP's Cloudflare challenge; that challenge
  appears probabilistic — a failed attempt is worth one immediate retry
  before concluding it's a hard block.
- **Taylor & Francis** (`tandfonline.com`): plain browser navigation
  worked without any Cloudflare gate at all, full text and PDF link both
  directly reachable.
- **Cowles Foundation Discussion Papers** (Phillips + coauthors): usually `cowles.yale.edu/sites/default/files/{YYYY-MM}/d{number}.pdf`, but the date folder isn't guessable from the DP number alone — two failed guesses this pass (`d2110`, and an early wrong folder for others) before finding the right one via the DP's `ideas.repec.org/p/cwl/cwldpp/{number}.html` page, which lists the exact working URL. Check that page first rather than guessing the folder.
- **Phillips's own homepage**: `korora.econ.yale.edu/phillips/pubs/art/p{number}.pdf` — direct PDFs, no landing page needed, but the server refuses HTTPS entirely (plain HTTP works).
- **Essex Research Repository** (Harvey/Leybourne/Taylor/Zu circle): `repository.essex.ac.uk` — CC-BY items download directly even when the landing page itself won't load for other tools.
- **EconStor** (`econstor.eu`, German National Library repository): reliable for German-affiliated authors (Hafner, Breitung) — find the exact bitstream path via the `econpapers.repec.org/RePEc:zbw:...` record page's redirect link, don't guess the filename.
- **MPRA** (`mpra.ub.uni-muenchen.de`): direct pattern `/{id}/1/MPRA_paper_{id}.pdf` worked as guessed.
- **Nottingham Granger Centre** (old-style hosting, `nottingham.ac.uk/research/groups/grangercentre/documents/`): works directly. Contrast with **Nottingham Repository** (`nottingham-repository.worktribe.com`, the newer institutional-repository system) — Cloudflare-gated with no Wayback fallback, see above. Same author circle, different hosting, very different accessibility.
- **SSRN**: blocks essentially all automated access now, even the plain abstract page (403) — not just the PDF. Don't spend time on SSRN URLs; go straight to looking for an institutional-repository or EconStor/MPRA mirror instead.
- **arXiv** (`econ.EM` category): reliably open, but `pdftotext` mangles subscript/superscript-heavy formulas often enough that a PyMuPDF page-render is the standing verification step before shipping any formula from any source, arXiv included.
- **JTSA 46(5) table of contents**: not reproducible from the paywalled Wiley ToC page directly (402 on fetch) — pulled instead via `ideas.repec.org/s/bla/jtsera.html`. Individual JTSA 46(5) articles, though, often have an open working-paper twin (Tinbergen, Lancaster, EconStor, Essex) findable by searching the paper title directly rather than going through Wiley at all.
