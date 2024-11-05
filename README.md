# surveillance: Temporal and Spatio-Temporal Modeling and Monitoring of Epidemic Phenomena

The open-source [R](https://www.R-project.org/) package **surveillance**
implements statistical methods for the *modeling* and *monitoring* of
epidemic phenomena based on (infectious disease) surveillance data.
This includes time series of counts, proportions and
categorical data as well as spatio-temporal point processes.
Potential users are biostatisticians, epidemiologists and others working in,
e.g., applied infectious disease epidemiology. However, applications could just
as well originate from environmetrics, reliability engineering, econometrics or
the social sciences.


## Prospective outbreak detection

[Salmon et al. (2016)](https://doi.org/10.18637/jss.v070.i10)
provide an overall guide to the monitoring capabilities of **surveillance**.
The paper is available as
[`vignette("monitoringCounts")`](https://CRAN.R-project.org/package=surveillance/vignettes/monitoringCounts.pdf)
with the package.
Further descriptions can be found in a book chapter by
H&ouml;hle and Mazick (2010,
[preprint](https://staff.math.su.se/hoehle/pubs/hoehle_mazick2009-preprint.pdf)),
and -- slightly outdated --
[H&ouml;hle (2007)](https://doi.org/10.1007/s00180-007-0074-8) or
[`vignette("surveillance")`](https://CRAN.R-project.org/package=surveillance/vignettes/surveillance.pdf).

- Aberration detection in *count data time series*,
  e.g., `farringtonFlexible()`.

- Online change-point detection in *categorical time series*,
  e.g., `categoricalCUSUM()`.\
  A Markov Chain approximation for computing the run-length distribution of the
  proposed likelihood ratio CUSUMs is available as function `LRCUSUM.runlength()`.

See the
[online reference index](https://surveillance.R-forge.R-project.org/pkgdown/reference/index.html#prospective-outbreak-detection)
for the complete list of algorithms.


## Modeling reporting delays

- Backprojection methods: `backprojNP()`

- Adjusting for occurred-but-not-yet-reported events: `nowcast()`, `bodaDelay()`

<!--
These methods have also been used for COVID-19 surveillance, for example:

- [*Analysis of the early COVID-19 epidemic curve in Germany by regression models
  with change points*](https://doi.org/10.1017/S0950268821000558)

- [*Nowcasting the COVID-19 pandemic in Bavaria*](https://doi.org/10.1002/bimj.202000112)
-->


## Endemic-epidemic modeling

[Meyer et al. (2017)](https://doi.org/10.18637/jss.v077.i11) provide a
guide to the spatio-temporal modeling capabilities of **surveillance**.
These so-called *endemic-epidemic* models have proven useful in a wide
range of applications, also beyond epidemiology.
A list of corresponding publications is maintained at
<https://surveillance.R-forge.R-project.org/applications_EE.html>.

### `twinstim()`
- models a *spatio-temporal point pattern* of infective events
- is described in
  [`vignette("twinstim")`](https://CRAN.R-project.org/package=surveillance/vignettes/twinstim.pdf)
- needs data of class `"epidataCS"`, which
  holds the observed events (with covariates) and exogenous covariates on a
  space-time grid (for the endemic/background component)
- features a model-based `epitest()` for space-time interaction
    
### `twinSIR()`
- models the susceptible-infectious-recovered
  *(SIR) event history of a fixed population*
- is described in
  [`vignette("twinSIR")`](https://CRAN.R-project.org/package=surveillance/vignettes/twinSIR.pdf)
- needs data of class `"epidata"`

### `hhh4()`
- models a (multivariate) *time series of infectious disease counts*
- is described in
  [`vignette("hhh4_spacetime")`](https://CRAN.R-project.org/package=surveillance/vignettes/hhh4_spacetime.pdf)
  for areal time series, and more generally in
  [`vignette("hhh4")`](https://CRAN.R-project.org/package=surveillance/vignettes/hhh4.pdf),
  including the univariate case
- needs data of class `"sts"` (see below)


## Data class `"sts"`

The S4 class `"sts"` (surveillance time series), created via `sts()` or
`linelist2sts()`, represents (multivariate) time series of counts.
For areal time series, the class can also capture population fractions, a map,
and a weight matrix.

For evaluation purposes, the package contains several datasets derived from the
[SurvStat@RKI](https://survstat.rki.de/) database maintained by the Robert Koch
Institute in Germany.
See the [online reference index](https://surveillance.R-forge.R-project.org/pkgdown/reference/index.html#datasets)
for the complete list of datasets.


## Installation

The stable release version of **surveillance** is hosted on the
Comprehensive R Archive Network (CRAN) at
<https://CRAN.R-project.org/package=surveillance>
and can be installed via

```R
install.packages("surveillance")
```

The development version of **surveillance** is hosted on R-Forge at
<https://R-Forge.R-project.org/projects/surveillance/> in a Subversion
(SVN) repository. It can be installed via

```R
install.packages("surveillance", repos = "https://R-Forge.R-project.org")
```

Alternatively, a development build can be installed from the
[R-universe](https://r-forge.r-universe.dev/surveillance) mirror of R-Forge.


## Feedback

Contributions are welcome!
Please report bugs via e-mail to `maintainer("surveillance")`.

Note that (large) new features are unlikely to be included in **surveillance**.
Some extensions have already been developed in separate packages, for example
[**hhh4contacts**](https://CRAN.R-project.org/package=hhh4contacts),
[**HIDDA.forecasting**](https://hidda.github.io/forecasting/),
[**hhh4addon**](https://github.com/jbracher/hhh4addon),
and [**hhh4ZI**](https://github.com/Junyi-L/hhh4ZI).


## Funding

The authors acknowledge financial support from the following institutions:

- [German Research Foundation](https://www.dfg.de/en/)
  (DFG, 2024--2027, #528691398)
- [FAU Interdisciplinary Center for Clinical Research](https://www.izkf.med.fau.de/en/)
  (IZKF, 2018--2021, junior project J75)
- [Swedish Research Council](https://www.vr.se/english.html)
  (VR, 2016--2019, #2015-05182)
- [Swiss National Science Foundation](https://www.snf.ch/en/)
  (SNSF, 2007--2010 and 2012--2015, projects
  [#116776](https://data.snf.ch/grants/grant/116776),
  [#124429](https://data.snf.ch/grants/grant/124429), and
  [#137919](https://data.snf.ch/grants/grant/137919))
- [Munich Center of Health Sciences](https://www.en.mc-health.uni-muenchen.de/)
  (MC-Health, 2007--2010)
- [German Research Foundation](https://www.dfg.de/en/) (DFG, 2003--2006)


## License

The **surveillance** package is free and open-source software,
and you are welcome to redistribute it under the terms of the
[GNU General Public License, version 2](https://www.gnu.org/licenses/gpl-2.0.html).
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY.
