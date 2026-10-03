---
title: "References, full bibliography"
blurb: "Full bibliography behind the replication record, organised by methodological family."
order: 7
---
This file lists every paper the project has touched, organized to match the
family files (`volatility-robustness.md`, `dating-and-root-inference.md` and
so on, see `README.md`). Those files hold the per-item detail, the formulas
and the exact-number verification, and this file is the index. The "Access"
column records what we actually found during research, and says nothing
general about a paper's open-access status.

the paper library holds local PDF copies. As of 2026-08-09, all 47 papers have a
local copy, linked inline. Thirty-one came through open-access routes and 14
through institutional access, once the session was routed through a UK
academic network (Jisc/Lancaster). "What changed between the first and second
pass" below explains the mechanics for each publisher. The last two, HLS 2017
and HLW 2020 (both Elsevier), resisted every automated route, including
Playwright on the institutional network, as "The two that automation couldn't
get" below describes. On 2026-08-11 we added four more specifically for the
DGP survey in [simulation-dgps.md](/replication/simulation-dgps), all open access, which
brings the total to 52.

## Foundational (pre-existing in exuber)

| Citation | Access | Role |
|---|---|---|
| Phillips, P.C.B., Wu, Y. & Yu, J. (2011). "Explosive Behavior in the 1990s Nasdaq: When Did Exuberance Escalate Asset Values?" *International Economic Review*, 52(1), 201-226. Working paper: Cowles Foundation DP 1699. | open | PWY sup-ADF, `radf()`'s `sadf`/`badf` |
| Phillips, P.C.B., Shi, S. & Yu, J. (2015a). "Testing for Multiple Bubbles: Historical Episodes of Exuberance and Collapse in the S&P 500." *International Economic Review*, 56(4), 1043-1078. | open | PSY/GSADF, `radf()`'s `gsadf` |
| Phillips, P.C.B., Shi, S. & Yu, J. (2015b). "Testing for Multiple Bubbles: Limit Theory of Real-Time Detectors." *International Economic Review*, 56(4), 1079-1134. Working paper: Cowles Foundation DP 1915. | open | BSADF sequence + `datestamp()` |
| Pavlidis et al. (2016) (cited by exuber's existing code, not re-verified) | not chased | sieve bootstrap already in `radf_sb.R`/`radf_sb_cv()` |

## Volatility robustness

| # | Citation | Access | Status |
|---|---|---|---|
| 1 | Harvey, D.I., Leybourne, S.J. & Zu, Y. (2020). "Sign-based unit root tests for explosive financial bubbles in the presence of deterministically time-varying volatility." *Econometric Theory*, 36(1), 122-169. `doi:10.1017/S0266466619000057` | open (institutional access) via Cambridge Core. Its 2025 OBES level-shift extension also recovered | todo |
| 2 | Kurozumi, E., Skrobotov, A. & Tsarev, A. (2024). "Time-Transformed Test for Bubbles under Non-stationary Volatility." *J. Financial Econometrics*. `doi:10.1093/jjfinec/nbae026` | open (arXiv:2012.13937) | done |
| 3 | Harvey, D.I., Leybourne, S.J. & Zu, Y. (2019). "Testing explosive bubbles with time-varying volatility." *Econometric Reviews*, 38(10), 1131-1151. | open (Nottingham Granger Centre DP 18/05) | done |
| 4 | Harvey, D.I., Leybourne, S.J., Taylor, A.M.R. & Zu, Y. (2024/2025). "A new heteroskedasticity-robust test for explosive bubbles." *JTSA*, 46(5), 846-866. `doi:10.1111/jtsa.12784` | open (CC-BY) | done (with-intercept variant) |
| 5 | Hafner, C.M. (2020). "Testing for Bubbles in Cryptocurrencies with Time-Varying Volatility." *J. Financial Econometrics*, 18(2), 233-249. | open (EconStor/IRTG 1792 DP 2018-005, SSRN itself blocks automated fetch) | evaluated, not implemented |
| 6 | Pedersen, T.Q. & Montes Schütte, E.C. (2020). "Testing for Explosive Bubbles in the Presence of Autocorrelated Innovations." *J. Empirical Finance*, 58, 207-225. | open (Aarhus CREATES RP 2017-9) | evaluated, mostly already covered by `radf_sb_cv()` |
| 7 | Kurozumi, E. & Nishi, M. (2025). "Testing for a bubble with a stochastically varying explosive coefficient." *JTSA*, 46(5), 945-965. | open (Open Access article); the article was never paywalled, but Wiley needs JS to serve it, which plain `curl` cannot do | todo, newly surfaced |
| 8 | Sarkar, A. & Wells, M.T. (2026). "Is There an AI Bubble? Robust Date-Stamping for Periods of Exuberance." arXiv:2604.12062 (preprint). Theory: Sarkar & Wells, "Double Local-to-Unity," arXiv:2512.06823. | open + | todo, preprint |
| 9 | Monschang, V. & Wilfling, B. (2021). "Sup-ADF-style bubble-detection methods under test." *Empirical Economics*, 61, 145-172. `doi:10.1007/s00181-020-01859-7` | open (CQE WP 78/2019) | new (2026-08-11), a TGARCH leverage and lognormal-mixture bubble DGP, see [simulation-dgps.md](/replication/simulation-dgps#12-tgarch11-with-leverage-effect) |
| 10 | Richter, S., Wang, W. & Wu, W.B. (2023, orig. 2018). "A supreme test for periodic explosive GARCH." *Econometrics* (MDPI); arXiv:1812.03475. | open | new (2026-08-11), a volatility-bubble GARCH (not a price-bubble one), out of scope for `radf()`, see [simulation-dgps.md](/replication/simulation-dgps#adjacent-not-a-price-level-dgp) |

## Dating and root inference

| Citation | Access | Notes |
|---|---|---|
| Harvey, D.I., Leybourne, S.J. & Sollis, R. (2017). "Improving the accuracy of asset price bubble start and end date estimators." *J. Empirical Finance*, 40, 121-138. `doi:10.1016/j.jempfin.2016.11.001` | open (supplied by the user, because automated fetching never got past the Elsevier/Nottingham Cloudflare gate; see the roundup below) | HLS |
| Harvey, D.I., Leybourne, S.J. & Whitehouse, E.J. (2020). "Date-stamping multiple bubble regimes." *J. Empirical Finance*, 58, 226-246. `doi:10.1016/j.jempfin.2020.06.004` | open (supplied by user; same blocker as HLS) | HLW |
| Pang, T., Du, L. & Chong, T.T.L. (2021). "Estimating multiple breaks in nonstationary autoregressive models." *J. Econometrics*, 221(1), 277-311. | open (MPRA 92074) | PDC |
| Kurozumi, E. & Skrobotov, A. (2023). "On the asymptotic behavior of bubble date estimators." *JTSA*, 44(4), 359-373. | open (arXiv:2110.04500) | KS |
| Kurozumi, E. & Skrobotov, A. (2023). "Improving the accuracy of bubble date estimators under time-varying volatility." arXiv:2306.02977. | open | newly found; a two-step WLS-based dating estimator |
| Kejriwal, M., Nguyen, L. & Perron, P. (2025). "An Improved Procedure for Retrospectively Dating the Emergence and Collapse of Bubbles." *JTSA*, 46(5). `doi:10.1111/jtsa.12810` | open (institutional access) | newly surfaced |
| Kurozumi, E. & Skrobotov, A. (2025). "Confidence Sets for the Emergence, Collapse, and Recovery Dates of a Bubble." arXiv:2511.16172 | open | newly surfaced |
| Phillips, P. C. B., & Magdalinos, T. (2007). "Limit theory for moderate deviations from a unit root." *J. Econometrics*, 136(1), 115-130. Working paper: Cowles Foundation DP 1471 (July 2004). | open Theorem 4.3 verified by PNG page-render 2026-08-09 | source-verified, not yet coded |
| Guo, G., Sun, Y. & Wang, S. (2019). "Testing for moderate explosiveness." *The Econometrics Journal*, 22(3), 279-303. | open (institutional access); we implemented `rootstamp()` from Skrobotov's secondary restatement before finding this paper, so the implementation is worth cross-checking against the primary source now that it is available | done, `rootstamp()` |
| Phillips, P.C.B., Magdalinos, T. & Giraitis, L. (2010). "Smoothing local-to-moderate unit root theory." *J. Econometrics*, 158(2), 274-279. Working paper: Cowles Foundation DP 1659. | open | background theory only |
| Phillips, P.C.B. & Shi, S. (2014). "Financial Bubble Implosion." Working paper: Cowles Foundation DP 1967. Published as "Financial Bubble Implosion and Reverse Regression," *Econometric Theory*. | open | newly surfaced; reverse-BSADF recovery dating |

## Monitoring

| Citation | Access | Notes |
|---|---|---|
| Homm, U. & Breitung, J. (2012). "Testing for speculative bubbles in stock markets: a comparison of alternative methods." *J. Financial Econometrics*, 10(1), 198-231. `doi:10.1093/jjfinec/nbr009` | open (institutional access) | formulas now directly readable, and not yet re-verified against this copy |
| Astill, S., Harvey, D.I., Leybourne, S.J., Taylor, A.M.R. & Zu, Y. (2021/2023). "CUSUM-Based Monitoring for Explosive Episodes in Financial Data in the Presence of Time-Varying Volatility." *J. Financial Econometrics*, 21(1), 187-227. `doi:10.1093/jjfinec/nbab009` | open (institutional access) | AHLTZ; the primary source is now available and we have not yet re-verified against this copy |
| Astill, S., Harvey, D.I., Leybourne, S.J., Sollis, R. & Taylor, A.M.R. (2018). "Real-Time Monitoring for Explosive Financial Bubbles." *JTSA*, 39, 863-891. | open (institutional access) | AHLST |
| Whitehouse, E.J., Harvey, D.I. & Leybourne, S.J. (2025). "Real-time monitoring procedures for early detection of bubbles." *Intl J. Forecasting*, 41(3), 1260-1277. `doi:10.1016/j.ijforecast.2024.12.005` | open (CC-BY) | supplied the verified AHLST decision rule and the FPR formula |
| Kurozumi, E. (2020). "Asymptotic properties of bubble monitoring tests." *Econometric Reviews*, 39(5), 510-538. `doi:10.1080/07474938.2019.1697086` | open (institutional access) | now readable at the primary source, where before we had only the abstract |
| Kurozumi, E. (2021). "Asymptotic Behavior of Delay Times of Bubble Monitoring Tests." *JTSA*, 42(3), 314-337. `doi:10.1111/jtsa.12569` | open (institutional access) | now readable at the primary source, where before we had only the abstract |
| Horváth, L. & Trapani, L. (2026). "Real-time monitoring with RCA models." *Econometric Theory*, 42, 514-547. Working paper arXiv:2312.11710. | open | |
| Astill, S., Taylor, A.M.R. & Zu, Y. (2026, forthcoming). "Covariate Augmented CUSUM Bubble Monitoring Procedures." *Econometric Theory*. Working paper: Essex Finance Centre WP 94. | open | the most useful single document we found for this file |
| Breitung, J. & Diegel, M. (2025). "Sequential Detector Statistics for Speculative Bubbles." *JTSA*, 46(5). | open (via EconStor) | newly surfaced, and a direct match |

## Multivariate

| Citation | Access | Notes |
|---|---|---|
| Chen, Y., Phillips, P.C.B. & Shi, S. (2020/2023). "Common Bubble Detection in Large Dimensional Financial Systems." *J. Financial Econometrics*, 21(4), 989-1063. Working paper: Cowles Foundation DP 2251. | open | done, `radf_common()` |
| Evripidou, A.C., Harvey, D.I., Leybourne, S.J. & Sollis, R. (2022). "Testing for Co-explosive Behaviour in Financial Time Series." *OBES*, 84(3), 624-650. `doi:10.1111/obes.12487` | open (institutional access) | evaluated; it was blocked and is now unblocked, so it is worth re-reading from the primary source |
| Phillips, P.C.B. & Yu, J. (2011). "Dating the Timeline of Financial Bubbles During the Subprime Crisis." *Quantitative Economics*, 2(3), 455-491. `doi:10.3982/QE82` | open (Cowles DP 1770) | evaluated, no new code needed |
| Greenaway-McGrevy, R. & Phillips, P.C.B. (2015/2016). "Hot Property in New Zealand: Empirical Evidence of Housing Bubbles in the Metropolitan Centres." *NZ Economic Papers*, 50(1), 88-113. Working paper: Cowles Foundation DP 2004. | open | evaluated, not implemented |

## Open research directions

| Citation | Access | Notes |
|---|---|---|
| Skrobotov, A. (2023). "Testing for explosive bubbles: a review." *Dependence Modeling*, 11(1), 1-26. | open (arXiv:2207.08249) | the source of the whole project |
| Lui, Y.L., Phillips, P.C.B. & Yu, J. (2024). "Robust Testing for Explosive Behavior with Strongly Dependent Errors." *J. Econometrics*, 238(2), 105626. Working paper: Cowles Foundation DP 2350. | open | evaluated, not implemented |
| Qian, J. & Su, L. (2016). "Shrinkage Estimation of Regression Models With Multiple Structural Changes." *Econometric Theory*, 32(6), 1376-1433. `doi:10.1017/S0266466615000237` | open (institutional access) via Cambridge Core; not itself a bubble-specific application | no bubble-specific LASSO precedent found |

## Alternative paradigms

| Citation | Access | Notes |
|---|---|---|
| Pavlidis, E.G. (2025). "Bubbles and crashes: A tale of quantiles." *JTSA*, 46(5), 884-907. | open (Lancaster ePrints) | quantile-based |
| Wu, R., Shi, S. & Wu, J. (2025). "Quantile analysis for financial bubble detection and surveillance." *JTSA*, 46(5), 908-931. | open (institutional access); SSRN itself is still blocked, so this copy came via Wiley instead | quantile-based |
| Blasques, F., Koopman, S.J., Mingoli, G. & Telg, S. (2025). "A Novel Test for the Presence of Local Explosive Dynamics." *JTSA*, 46(5), 966-980. `doi:10.1111/jtsa.70001` | open (Tinbergen DP 24-036/III) | noncausal |
| Gourieroux, C. & Jasiak, J. (2025). "A Stochastic Tree for Bubble Asset Modelling and Pricing." *JTSA*, 46(5), 932-944. | open (institutional access, and the article is Open Access itself) | out of scope, since it concerns pricing and not testing |
| Bhandari, A. "Rational Bubbles at the Spectral Edge." arXiv:2607.03933. | open | out of scope, since it is an entirely different paradigm |
| Chan, J.C.C. & Santi, C. (2021). "Speculative Bubbles in Present-Value Models: A Bayesian Markov-Switching State Space Approach." *J. Economic Dynamics and Control*, 127, 104101. | open (author's site) | new (2026-08-11), Markov-switching bubble/collapse timing, see [simulation-dgps.md](/replication/simulation-dgps#14-markov-switching-present-value-bubble) |
| Chen, H., Chen, L., Huang, D., Li, Y. & Zhang, Z. (2026). "Technology Fundamentals and False Bubble Detection: Evidence from Dot-Com and AI Episodes." arXiv:2604.25826. | open | new (2026-08-11), deterministic false-bubble null, see [simulation-dgps.md](/replication/simulation-dgps#15-deterministic-technology-adoption-false-bubble-null) |

## Practitioner guidance

| Citation | Access |
|---|---|
| Phillips, P.C.B., Shi, S. & Yu, J. (2014). "Specification Sensitivity in Right-Tailed Unit Root Testing for Explosive Behaviour." *OBES*, 76(3), 315-333. | open |
| Phillips, P.C.B. & Shi, S. (2019). "Detecting Financial Collapse and Ballooning Sovereign Risk." *OBES*, 81(6), 1336-1361. Working paper: Cowles Foundation DP 2110. | open (institutional access); SSRN and the SMU repository copy are both still blocked, so this copy came via Wiley instead |
| Basele, R.B., Phillips, P.C.B. & Shi, S. (2025). "Speculative Bubbles in the Recent AI Boom: Nasdaq and the Magnificent Seven." *JTSA*, 46(5). | open |

## Context / meta

- Harvey, D.I. & Leybourne, S.J. (2025). Guest editors' introduction, *JTSA* special issue "Recent Developments in Time-Series Methods for Detecting Bubbles and Crashes," 46(5). `doi:10.1111/jtsa.70003`. The editorial itself is paywalled, but its existence is the reason so much of this project received a 2025-26 refresh. We recovered the full issue table of contents through RePEc/IDEAS and not the paywalled editorial text.
- Hu, Y. et al. (2023). "A review of Phillips-type right-tailed unit root bubble detection tests." *J. Economic Surveys*, 37(1), 141-158. This is a second review alongside Skrobotov (2023). We have not yet used it as a source for this project, and it is a candidate for a second pass.

## The two that automation couldn't get: HLS and HLW

Harvey, Leybourne & Sollis (2017) and Harvey, Leybourne & Whitehouse (2020)
(both *J. Empirical Finance*, in
[dating-and-root-inference.md](/replication/dating-and-root-inference)) were the last
two holdouts. Both are in the paper library anyway, because the user supplied them
directly and no tool here fetched them. It is worth recording why the
automated routes failed, since this is the one clean counterexample in an
otherwise near-total recovery. Both papers are Elsevier (ScienceDirect), and
Unpaywall lists an open Green-OA copy of each at
`nottingham-repository.worktribe.com`. That repository returns 403 to every
request behind a Cloudflare "Just a moment…" challenge, and the challenge held
even against a real Playwright browser session routed through the same
academic network that got everything else through. The ScienceDirect article
pages hit the identical wall. The old-style Nottingham Granger Centre hosting
serves files directly, with no gate at all. The newer institutional-repository
platform evidently runs a stricter bot check than Wiley, OUP, Taylor & Francis
or Cambridge Core did, and none of those blocked the same browser session. The
Wayback Machine has only ever crawled the HTML landing pages for these two
papers and never the underlying PDF asset (we checked via the CDX API), so that
fallback did not apply either.

## What changed between the first and second pass

The first pass (plain `curl`, no institutional network) found 31 of 47 papers
and left 16 "not found". The second pass ran once the session was routed
through the Jisc/Lancaster academic network. It recovered 14 of those 16 by
the same DOI, as follows.

- Cambridge Core worked immediately with plain `curl` once it recognized the
  IP address. It needed no JS and no browser, and 2 of 2 papers were
  recovered this way.
- Wiley (`onlinelibrary.wiley.com`) is a JS-rendered SPA regardless of
  subscription status, and `curl` gets a 60KB Angular shell whatever we do.
  We needed a real browser (Playwright) to load the page. We then called
  `page.evaluate(() => fetch(...))` inside that page's JS context to pull the
  PDF bytes as base64, because Chrome's built-in PDF viewer intercepts a plain
  navigation before the raw bytes can be captured any other way. Eight of
  eight papers were recovered this way. One of them, Kurozumi & Nishi 2025,
  turned out to be Open Access and was never paywalled.
- Oxford University Press (`academic.oup.com`) and Taylor & Francis
  (`tandfonline.com`) also needed a real browser. The first attempts hit the
  same Cloudflare challenge as Elsevier and Nottingham, but later attempts
  against the same URLs got through cleanly, which suggests a probabilistic or
  rate-limited challenge and no hard block. The PDFs were then on a different
  subdomain (`watermark02.silverchair.com`, the PDF host of OUP), and the same
  `fetch()`-from-page-context trick worked there too. Four of four papers were
  recovered this way, after retries.
- Elsevier (`sciencedirect.com`) and the Nottingham repository never got past
  the Cloudflare challenge, with or without a browser, as described above.

The net effect is less "here are more papers" than "many formulas that we had
sourced secondhand (from restatements, abstracts or Skrobotov's review) can
now be verified against the actual primary text". Several rows above say
"worth re-reading from primary source now that it's available", and that is
real follow-up work and more than a better citation list.

## Recurring access routes (useful for future passes)

- On an institutional network, prefer the publisher's own page over hunting
  for a working-paper mirror. Once the network was recognized, Wiley, OUP,
  T&F and Cambridge Core all served full text directly, with none of the
  Cowles/arXiv/MPRA/EconStor detective work described below. That work
  matters mainly when we are not on such a network, and for the few
  publishers (Elsevier, and repositories fronted by the newer Cloudflare-grade
  bot checks) that block automated access whatever the subscription status.
- Wiley (`onlinelibrary.wiley.com`) needs a real browser, because it is a
  JS-rendered SPA for everyone, subscriber or not. The endpoint is
  `doi/pdfdirect/{doi}`. Fetch it with `page.evaluate(() => fetch(...))` from
  within an already-loaded Wiley page (same origin, so no CORS issue),
  base64-encode the bytes in the page, and decode them outside the browser.
  Navigating to the URL directly instead hands the bytes to the built-in PDF
  viewer extension of Chrome, which intercepts them before they can be
  captured as a plain network response body.
- OUP (`academic.oup.com`) and its PDF host Silverchair
  (`watermark02.silverchair.com`) accept the same `fetch()`-from-page-context
  technique once we are past the Cloudflare challenge of OUP. That challenge
  appears probabilistic, so a failed attempt deserves one immediate retry
  before we conclude it is a hard block.
- Taylor & Francis (`tandfonline.com`) worked with plain browser navigation
  and no Cloudflare gate at all. Both the full text and the PDF link were
  directly reachable.
- Cowles Foundation Discussion Papers (Phillips and coauthors) usually live at
  `cowles.yale.edu/sites/default/files/{YYYY-MM}/d{number}.pdf`, but the date
  folder cannot be guessed from the DP number alone. We made two failed guesses
  in this pass (`d2110`, and an early wrong folder for others) before finding
  the right one on the page `ideas.repec.org/p/cwl/cwldpp/{number}.html` of the
  DP, which lists the exact working URL. Check that page first and do not guess
  the folder.
- Phillips's own homepage serves direct PDFs at
  `korora.econ.yale.edu/phillips/pubs/art/p{number}.pdf`, with no landing page
  needed. The server refuses HTTPS entirely, but plain HTTP works.
- The Essex Research Repository (`repository.essex.ac.uk`, home of the
  Harvey/Leybourne/Taylor/Zu circle) downloads CC-BY items directly, even when
  the landing page itself will not load for other tools.
- EconStor (`econstor.eu`, the repository of the German National Library) is
  reliable for German-affiliated authors (Hafner, Breitung). Find the exact
  bitstream path through the redirect link on the
  `econpapers.repec.org/RePEc:zbw:...` record page, and do not guess the
  filename.
- MPRA (`mpra.ub.uni-muenchen.de`) worked with the direct pattern
  `/{id}/1/MPRA_paper_{id}.pdf`, exactly as we guessed it.
- Nottingham Granger Centre (old-style hosting,
  `nottingham.ac.uk/research/groups/grangercentre/documents/`) works
  directly. The Nottingham Repository (`nottingham-repository.worktribe.com`,
  the newer institutional-repository system) is Cloudflare-gated with no
  Wayback fallback, as described above. The author circle is the same, but the
  hosting and the accessibility differ completely.
- SSRN now blocks essentially all automated access, including the plain
  abstract page (403) and not only the PDF. Do not spend time on SSRN URLs.
  Look instead for an institutional-repository, EconStor or MPRA mirror.
- arXiv (`econ.EM` category) is reliably open, but `pdftotext` often mangles
  formulas heavy in subscripts and superscripts. A PyMuPDF page render is
  therefore the standing verification step before we ship any formula from any
  source, arXiv included.
- The table of contents of JTSA 46(5) cannot be reproduced from the paywalled
  Wiley ToC page directly (402 on fetch), so we pulled it from
  `ideas.repec.org/s/bla/jtsera.html`. Individual JTSA 46(5) articles often
  have an open working-paper twin (Tinbergen, Lancaster, EconStor, Essex),
  which we can find by searching for the paper title and without going through
  Wiley at all.
