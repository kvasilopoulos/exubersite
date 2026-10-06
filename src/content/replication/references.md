---
title: "References, full bibliography"
blurb: "Full bibliography behind the replication record, organised by methodological family."
order: 8
---
This file lists every paper the project has touched, organized to match the
family files (`volatility-robustness.md`, `dating-and-root-inference.md` and
so on, see `README.md`). Those files hold the per-item detail, the formulas
and the exact-number verification, and this file is the index.

The "Access" column says where a copy of each paper can be obtained.

the paper library holds local PDF copies of all 52 papers, linked inline. Harvey, Leybourne & Sollis (2017) and Harvey, Leybourne & Whitehouse (2020), both in the *Journal of Empirical Finance*, are Elsevier papers and are not freely available from the publisher.

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
| 1 | Harvey, D.I., Leybourne, S.J. & Zu, Y. (2020). "Sign-based unit root tests for explosive financial bubbles in the presence of deterministically time-varying volatility." *Econometric Theory*, 36(1), 122-169. `doi:10.1017/S0266466619000057` | open (institutional access) via Cambridge Core. Its 2025 OBES level-shift extension also recovered | |
| 2 | Kurozumi, E., Skrobotov, A. & Tsarev, A. (2024). "Time-Transformed Test for Bubbles under Non-stationary Volatility." *J. Financial Econometrics*. `doi:10.1093/jjfinec/nbae026` | open (arXiv:2012.13937) | done |
| 3 | Harvey, D.I., Leybourne, S.J. & Zu, Y. (2019). "Testing explosive bubbles with time-varying volatility." *Econometric Reviews*, 38(10), 1131-1151. | open (Nottingham Granger Centre DP 18/05) | done |
| 4 | Harvey, D.I., Leybourne, S.J., Taylor, A.M.R. & Zu, Y. (2024/2025). "A new heteroskedasticity-robust test for explosive bubbles." *JTSA*, 46(5), 846-866. `doi:10.1111/jtsa.12784` | open (CC-BY) | done (with-intercept variant) |
| 5 | Hafner, C.M. (2020). "Testing for Bubbles in Cryptocurrencies with Time-Varying Volatility." *J. Financial Econometrics*, 18(2), 233-249. | open (EconStor/IRTG 1792 DP 2018-005, SSRN itself blocks automated fetch) | evaluated, not implemented |
| 6 | Pedersen, T.Q. & Montes Schütte, E.C. (2020). "Testing for Explosive Bubbles in the Presence of Autocorrelated Innovations." *J. Empirical Finance*, 58, 207-225. | open (Aarhus CREATES RP 2017-9) | evaluated, mostly already covered by `radf_sb_cv()` |
| 7 | Kurozumi, E. & Nishi, M. (2025). "Testing for a bubble with a stochastically varying explosive coefficient." *JTSA*, 46(5), 945-965. | open (Open Access article); the article was never paywalled, but Wiley needs JS to serve it, which plain `curl` cannot do | todo, newly surfaced |
| 8 | Sarkar, A. & Wells, M.T. (2026). "Is There an AI Bubble? Robust Date-Stamping for Periods of Exuberance." arXiv:2604.12062 (preprint). Theory: Sarkar & Wells, "Double Local-to-Unity," arXiv:2512.06823. | open + | |
| 9 | Monschang, V. & Wilfling, B. (2021). "Sup-ADF-style bubble-detection methods under test." *Empirical Economics*, 61, 145-172. `doi:10.1007/s00181-020-01859-7` | open (CQE WP 78/2019) | DGP source, a TGARCH leverage and lognormal-mixture bubble DGP, see [simulation-dgps.md](/replication/simulation-dgps#12-tgarch11-with-leverage-effect) |
| 10 | Richter, S., Wang, W. & Wu, W.B. (2023, orig. 2018). "A supreme test for periodic explosive GARCH." *Econometrics* (MDPI); arXiv:1812.03475. | open | DGP source, a volatility-bubble GARCH (not a price-bubble one), out of scope for `radf()`, see [simulation-dgps.md](/replication/simulation-dgps#adjacent-not-a-price-level-dgp) |

## Dating and root inference

| Citation | Access | Notes |
|---|---|---|
| Harvey, D.I., Leybourne, S.J. & Sollis, R. (2017). "Improving the accuracy of asset price bubble start and end date estimators." *J. Empirical Finance*, 40, 121-138. `doi:10.1016/j.jempfin.2016.11.001` | open (supplied by the user, because automated fetching never got past the Elsevier/Nottingham Cloudflare gate; see the roundup below) | HLS |
| Harvey, D.I., Leybourne, S.J. & Whitehouse, E.J. (2020). "Date-stamping multiple bubble regimes." *J. Empirical Finance*, 58, 226-246. `doi:10.1016/j.jempfin.2020.06.004` | open (supplied by user; same blocker as HLS) | HLW |
| Pang, T., Du, L. & Chong, T.T.L. (2021). "Estimating multiple breaks in nonstationary autoregressive models." *J. Econometrics*, 221(1), 277-311. | open (MPRA 92074) | PDC |
| Kurozumi, E. & Skrobotov, A. (2023). "On the asymptotic behavior of bubble date estimators." *JTSA*, 44(4), 359-373. | open (arXiv:2110.04500) | KS |
| Kurozumi, E. & Skrobotov, A. (2023). "Improving the accuracy of bubble date estimators under time-varying volatility." arXiv:2306.02977. | open | newly found; a two-step WLS-based dating estimator |
| Kejriwal, M., Nguyen, L. & Perron, P. (2025). "An Improved Procedure for Retrospectively Dating the Emergence and Collapse of Bubbles." *JTSA*, 46(5). `doi:10.1111/jtsa.12810` | open (institutional access) | |
| Kurozumi, E. & Skrobotov, A. (2025). "Confidence Sets for the Emergence, Collapse, and Recovery Dates of a Bubble." arXiv:2511.16172 | open | |
| Phillips, P. C. B., & Magdalinos, T. (2007). "Limit theory for moderate deviations from a unit root." *J. Econometrics*, 136(1), 115-130. Working paper: Cowles Foundation DP 1471 (July 2004). | open | Theorem 4.3 |
| Guo, G., Sun, Y. & Wang, S. (2019). "Testing for moderate explosiveness." *The Econometrics Journal*, 22(3), 279-303. | open (institutional access) | done, `rootstamp()` |
| Phillips, P.C.B., Magdalinos, T. & Giraitis, L. (2010). "Smoothing local-to-moderate unit root theory." *J. Econometrics*, 158(2), 274-279. Working paper: Cowles Foundation DP 1659. | open | background theory only |
| Phillips, P.C.B. & Shi, S. (2014). "Financial Bubble Implosion." Working paper: Cowles Foundation DP 1967. Published as "Financial Bubble Implosion and Reverse Regression," *Econometric Theory*. | open | newly surfaced; reverse-BSADF recovery dating |

## Monitoring

| Citation | Access | Notes |
|---|---|---|
| Homm, U. & Breitung, J. (2012). "Testing for speculative bubbles in stock markets: a comparison of alternative methods." *J. Financial Econometrics*, 10(1), 198-231. `doi:10.1093/jjfinec/nbr009` | open (institutional access) | |
| Astill, S., Harvey, D.I., Leybourne, S.J., Taylor, A.M.R. & Zu, Y. (2021/2023). "CUSUM-Based Monitoring for Explosive Episodes in Financial Data in the Presence of Time-Varying Volatility." *J. Financial Econometrics*, 21(1), 187-227. `doi:10.1093/jjfinec/nbab009` | open (institutional access) | AHLTZ |
| Astill, S., Harvey, D.I., Leybourne, S.J., Sollis, R. & Taylor, A.M.R. (2018). "Real-Time Monitoring for Explosive Financial Bubbles." *JTSA*, 39, 863-891. | open (institutional access) | AHLST |
| Whitehouse, E.J., Harvey, D.I. & Leybourne, S.J. (2025). "Real-time monitoring procedures for early detection of bubbles." *Intl J. Forecasting*, 41(3), 1260-1277. `doi:10.1016/j.ijforecast.2024.12.005` | open (CC-BY) | supplied the verified AHLST decision rule and the FPR formula |
| Kurozumi, E. (2020). "Asymptotic properties of bubble monitoring tests." *Econometric Reviews*, 39(5), 510-538. `doi:10.1080/07474938.2019.1697086` | open (institutional access) | |
| Kurozumi, E. (2021). "Asymptotic Behavior of Delay Times of Bubble Monitoring Tests." *JTSA*, 42(3), 314-337. `doi:10.1111/jtsa.12569` | open (institutional access) | |
| Horváth, L. & Trapani, L. (2026). "Real-time monitoring with RCA models." *Econometric Theory*, 42, 514-547. Working paper arXiv:2312.11710. | open | |
| Astill, S., Taylor, A.M.R. & Zu, Y. (2026, forthcoming). "Covariate Augmented CUSUM Bubble Monitoring Procedures." *Econometric Theory*. Working paper: Essex Finance Centre WP 94. | open | |
| Breitung, J. & Diegel, M. (2025). "Sequential Detector Statistics for Speculative Bubbles." *JTSA*, 46(5). | open (via EconStor) | newly surfaced, and a direct match |

## Multivariate

| Citation | Access | Notes |
|---|---|---|
| Chen, Y., Phillips, P.C.B. & Shi, S. (2020/2023). "Common Bubble Detection in Large Dimensional Financial Systems." *J. Financial Econometrics*, 21(4), 989-1063. Working paper: Cowles Foundation DP 2251. | open | done, `radf_common()` |
| Evripidou, A.C., Harvey, D.I., Leybourne, S.J. & Sollis, R. (2022). "Testing for Co-explosive Behaviour in Financial Time Series." *OBES*, 84(3), 624-650. `doi:10.1111/obes.12487` | open (institutional access) | evaluated |
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
| Wu, R., Shi, S. & Wu, J. (2025). "Quantile analysis for financial bubble detection and surveillance." *JTSA*, 46(5), 908-931. | open (institutional access) | quantile-based |
| Blasques, F., Koopman, S.J., Mingoli, G. & Telg, S. (2025). "A Novel Test for the Presence of Local Explosive Dynamics." *JTSA*, 46(5), 966-980. `doi:10.1111/jtsa.70001` | open (Tinbergen DP 24-036/III) | noncausal |
| Gourieroux, C. & Jasiak, J. (2025). "A Stochastic Tree for Bubble Asset Modelling and Pricing." *JTSA*, 46(5), 932-944. | open (institutional access, and the article is Open Access itself) | out of scope, since it concerns pricing and not testing |
| Bhandari, A. "Rational Bubbles at the Spectral Edge." arXiv:2607.03933. | open | out of scope, since it is an entirely different paradigm |
| Chan, J.C.C. & Santi, C. (2021). "Speculative Bubbles in Present-Value Models: A Bayesian Markov-Switching State Space Approach." *J. Economic Dynamics and Control*, 127, 104101. | open (author's site) | DGP source, Markov-switching bubble/collapse timing, see [simulation-dgps.md](/replication/simulation-dgps#14-markov-switching-present-value-bubble) |
| Chen, H., Chen, L., Huang, D., Li, Y. & Zhang, Z. (2026). "Technology Fundamentals and False Bubble Detection: Evidence from Dot-Com and AI Episodes." arXiv:2604.25826. | open | DGP source, deterministic false-bubble null, see [simulation-dgps.md](/replication/simulation-dgps#15-deterministic-technology-adoption-false-bubble-null) |

## Practitioner guidance

| Citation | Access |
|---|---|
| Phillips, P.C.B., Shi, S. & Yu, J. (2014). "Specification Sensitivity in Right-Tailed Unit Root Testing for Explosive Behaviour." *OBES*, 76(3), 315-333. | open |
| Phillips, P.C.B. & Shi, S. (2019). "Detecting Financial Collapse and Ballooning Sovereign Risk." *OBES*, 81(6), 1336-1361. Working paper: Cowles Foundation DP 2110. | open (institutional access) |
| Basele, R.B., Phillips, P.C.B. & Shi, S. (2025). "Speculative Bubbles in the Recent AI Boom: Nasdaq and the Magnificent Seven." *JTSA*, 46(5). | open |

## Context / meta

- Harvey, D.I. & Leybourne, S.J. (2025). Guest editors' introduction, *JTSA* special issue "Recent Developments in Time-Series Methods for Detecting Bubbles and Crashes," 46(5). `doi:10.1111/jtsa.70003`. The issue collects most of the recent work this project draws on.
- Hu, Y. et al. (2023). "A review of Phillips-type right-tailed unit root bubble detection tests." *J. Economic Surveys*, 37(1), 141-158. This is a second review alongside Skrobotov (2023). It is not used as a source in this project.
