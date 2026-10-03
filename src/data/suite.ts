// The three implementations and the same run in each language. The homepage's
// "get started" section and /guide (Getting started) both use this.

export const suite = [
  {
    name: "exubercore",
    lang: "C++ / Armadillo",
    status: "core" as const,
    description: "The recursive least-squares statistic itself. It has no R or Python dependency: it takes plain matrices and returns a statistic. Anything that needs random numbers (simulation, bootstrap, date-stamping) stays in the host language for now.",
    href: "https://github.com/kvasilopoulos/exubercore",
  },
  {
    name: "exuber",
    lang: "R · on CRAN",
    status: "stable" as const,
    description: "The complete package: the statistic, Monte Carlo and wild-bootstrap critical values, date-stamping, plotting and panel methods.",
    install: 'install.packages("exuber")',
    href: "https://github.com/kvasilopoulos/exuber",
  },
  {
    name: "pyexuber",
    lang: "Python · pybind11",
    status: "early" as const,
    description: "Bindings to the same core, imported as exuber. The R functions are ported method by method (tests, bootstraps, dating, monitoring, multivariate methods and simulation) and cross-checked against R.",
    install: "git clone …exuber/pyexuber && uv sync --dev",
    href: "https://github.com/kvasilopoulos/pyexuber",
  },
];

export const examples = [
  {
    label: "R",
    lang: "r",
    code: `library(exuber)

rsim <- radf(sim_data)               # adf, badf, sadf, gsadf, bsadf, ...
cv <- radf_mc_cv(n = NROW(sim_data)) # simulated 90/95/99% thresholds

summary(rsim, cv)   # rejects/does not reject H0 per series, per statistic
autoplot(rsim, cv)  # bsadf vs. its threshold, explosive spans shaded
datestamp(rsim, cv) # a Start/End/Duration row per detected episode`,
  },
  {
    label: "Python",
    lang: "python",
    code: `import numpy as np
import exuber

data = np.cumsum(np.random.randn(200))
result = exuber.radf(data)          # RadfResult: .adf .badf .sadf .gsadf .bsadf
cv = exuber.radf_mc_cv(n=len(data)) # RadfCv: simulated 90/95/99% thresholds

result.adf, result.sadf, result.gsadf
exuber.datestamp(result, cv)  # list[Episode]: one per detected explosive span`,
  },
  {
    label: "C++",
    lang: "cpp",
    code: `#include <exubercore/radf.hpp>

arma::vec y = arma::cumsum(arma::randn(200));
auto result = exubercore::radf(y, /* min_win = */ 20);
// arma::vec: badf, then adf, sadf, gsadf, then the bsadf sequence --
// see the reference entry for the exact layout and offsets.`,
  },
];
