# surveillance 1.25.0 (2025-06-24)

## New Features

- `algo.cusum()` gains a `reset` option: if enabled, the CUSUM statistic
  restarts from 0 after an alarm.  (Wish of Ann Christin Vietor.)

- `intensityplot.twinstim()` in its spatial variant now automatically
  labels the color key.

- `intensityplot.twinstim()` now also supports plotting the component
  intensities, not just their proportions.

- `plotHHH4_neweights()` now only excludes distance 0 by default when the
  model has an AR component.

- `plot.sts()` now allows `type = ~time` and `type = ~unit` as short forms
  of `type = observed ~ time` and `type = observed ~ unit`, respectively.

## Bug Fixes

- `algo.glrnb()` warns about unimplemented `dir`ection settings.

- `intensityplot.twinstim(..., aggregate = "space")` no longer disables
  `checkEmptyRC` when it calls `sp::spplot()`, so grids with
  horizontal or vertical gaps are now plotted without artifacts.

- `plot.hhh4()` can now be called with an *unnamed* `type` argument and
  additional named arguments of specific methods.

## Deprecated & Defunct

- The deprecated function `stsplot_spacetime()` has been removed.
  `plot.sts()` with the old `type = observed ~ 1 | unit` now uses
  `stsplot_space()` with a warning.


# surveillance 1.24.1 (2024-11-05)

This maintenance release adjusts a test of `formatDate()`
to be compatible with recent changes in R-devel.


# surveillance 1.24.0 (2024-10-01)

## New Features

- `residuals.hhh4()` can now also compute `type = "pearson"` residuals.

- The standard `frequency()` and `start()` generics can now be used to
  extract the respective slots from an `"sts"` object.

- When indexing a spatial `"sts"` object that contains a map by region names,
  i.e., `<sts>[,<character>]`, the additional argument `drop = TRUE` can
  now be used to subset the map accordingly.

## Package Infrastructure

- **sp** version 2.1-4 is now required, mainly to skip versions that
  produce misleading/obsolete startup messages or throw unnecessary
  warnings when **sf** is not available.

## Deprecated & Defunct

- The experimental `algo.twins()` implementation of
  [Held et al. (2006)](https://doi.org/10.1093/biostatistics/kxj016)
  has been removed from the package. The source code has been migrated
  to a separate R package, **twins**, which is archived at
  <https://codeberg.org/EE-hub/twins>.

## Bug Fixes

- `update.hhh4(S =)` can now be used for a previously disabled model
  component, when it simply uses `f = ~1` and adds the requested number
  of harmonics (but gives a warning as the specification may need to be
  tweaked, e.g., for an offset).


# surveillance 1.23.1 (2024-09-02)

## New Features

- The `sts()` constructor now also accepts an `"sf"` object as `map` input;
  it is internally converted to `"SpatialPolygons"` as required by the
  `"sts"` class.  (Based on a patch by Sophie Reichert.)

- `as.epidataCS()` is faster in determining potential event sources.

## Package Infrastructure

- **Rcpp** is no longer used. Only two small helper functions (for
  `backprojNP(eq3a.method="C")` and `as.epidataCS()`) were using it
  (inefficiently) and have been replaced by C implementations.
  This also reduces the size of the installed package.


# surveillance 1.23.0 (2024-05-03)

## New Features

- `update.epidataCS()` gained an argument `stgrid` to update the
  spatio-temporal grid data in an existing `"epidataCS"` object.
  This enables updates/transformations of endemic variables and/or changes
  of the time intervals without needing to do `as.epidataCS()` from scratch.

## Bug Fixes

- Start values for endemic intercepts in `twinstim()` are now robust
  against non-finite values in offset terms. 

- `intensityplot.twinstim(aggregate = "space")` no longer fails for
  endemic-only fits.


# surveillance 1.22.1 (2023-11-27)

## Bug Fixes

- The `pit()` plot could lack some tick marks on the y-axis (for R >= 4.2.0).


# surveillance 1.22.0 (2023-10-30)

## Package Infrastructure

- Legacy functions `unionSpatialPolygons()` and `polyAtBorder()` now use
  **sf** in place of **rgeos**.

## Deprecated & Defunct

- Long unused methods for `"gpc.poly"` objects (from package **gpclib**)
  have now been removed and `surveillance.options("gpclib")` is obsolete.

- Package **rgeos** is no longer available as a `clipper` method for
  `as.epidataCS()`. The previous default **polyclip** remains as the only
  option.


# surveillance 1.21.1 (2023-05-16)

This is a maintenance release, fixing encoding-related portability
issues and increasing test coverage for rarely used functionality.


# surveillance 1.21.0 (2023-03-14)

## Minor Changes

- `nbOrder()` has been re-implemented: it is now more efficient and no
  longer depends on **spdep**. Furthermore, it now defaults to
  `maxlag = Inf`; the historical default `maxlag = 1` was barely useful.
  It no longer messages (about the range of the detected orders).

- Printing `"sts"` objects with a map now shows the first row of the
  attached data (if present) instead of the object summary.

## Package Infrastructure

- Accommodate the current evolution of **sp**: **sf** is suggested and
  some examples are now conditionalized on its availability.

## Deprecated & Defunct

- **surveillance** no longer relies on the **maptools** package:
  `unionSpatialPolygons()` with `method = "gpclib"` is deprecated
  and now uses the default method with a warning.

- Long unused `scale.gpc.poly()` and `inside.gpc.poly()` are deprecated;
  the unused and undocumented `diameter.gpc.poly()` method has been removed.

- `stsplot_spacetime()` is formally deprecated; it has long been
  superseded by `stsplot_space()` and an `animate()` method for `"sts"`
  objects.


# surveillance 1.20.3 (2022-11-14)

## Package Infrastructure

- `vignette("monitoringCounts")` now uses **knitr** as its engine
  to work around [Bug 18318](https://bugs.R-project.org/show_bug.cgi?id=18318).


# surveillance 1.20.2 (2022-10-31)

## New Features

- `plotHHH4_fitted()` can now produce simple (unformatted) time indexes
  if argument `xaxis = NA`.

## Minor Changes

- Various documentation improvements,
  including an example for `predict.hhh4()`.

- `intensityplot.twinstim()` no longer depends on package
  [**maptools**](https://CRAN.R-project.org/package=maptools).

## Bug Fixes

- `hhh4()` now warns about interaction terms in model formulae.
  These are not implemented and were silently ignored previously.

- Fixed a memory leak in `algo.twins()`; note that this old experimental
  MCMC implementation for a two-component epidemic model may be removed in
  future versions of the package.


# surveillance 1.20.1 (2022-07-13)

## Bug Fixes

- `ks.plot.unif()`: accommodate to `NO_S_TYPEDEFS` in R >= 4.3.0.

- `boda()` with `samplingMethod="marginals"` gave all-`NA` upperbounds
  in INLA >= 21.07.10. `boda()` now also works around a scoping issue
  (with `E`) in recent versions of INLA that led to wrongly scaled
  upperbounds.


# surveillance 1.20.0 (2022-02-15)

## New Features

- `plotHHH4_season()` gained a `period` argument to support harmonics
  with periods longer than the frequency of the `"sts"` object.

- `stsplot_space()` now supports passing a `col` argument to `sp::spplot()`
  to change the colour of the polygon lines.

- `plotHHH4_fitted()` can now handle time series with missing values.

## Minor Changes

- If the Nelder-Mead optimizer is used for the variance parameters in
  `hhh4()`, it is now limited to 500 (not 300) iterations by default
  (consistent with the default in `optim()`).

- Printing an `"sts"` object now omits the `neighbourhood` component if
  that was not set (all-`NA` prototype).

- `simulate.hhh4(..., simplify = TRUE)` now consistently returns a
  3d array (nTime x nUnit x nsim), even for `nsim = 1` (for which plotting
  now works).

- The default legend in `stsplot_time1()` now only includes plotted elements.

- `wrap.algo()` no longer prints progress when there is only one area.

- `summary.hhh4()` now prints the number of excluded observations (due to
  missingness), if any.

## Bug Fixes

- The `print`-method for `summary.hhh4()` did not apply the `digits`
  argument to the coefficient matrix. Furthermore, printing of estimated
  variance parameters now adheres to *significant* `digits` as documented.

- The `[`-method for the `"hhh4sims"` class was not registered and thus
  only available internally. Array-like subsetting of simulated counts
  now retains the class.

- `farringtonFlexible()` with activated `populationOffset` (non-default)
  always used the population data of the *first* time series in the fitting step
  while iterating over a multivariate `"sts"` object.

- `plotHHH4_ri(..., exp = TRUE)` failed to use a log-scale color axis if
  further `colorkey` options were passed in a list. The (default) color breaks
  could fail to span the range of the data without warning (resulting in
  unfilled polygons). This is now checked and the default breaks are now
  equally spaced on the log-scale.

- `stsplot_time1()` did not pass `lty` to `polygon()`
  and `lwd` to `legend()`.

- `rps()` was wrong for distributions close to a point mass at zero,
  e.g., for `mu = 1e-3` and `x >= 4`. It is now also protected
  against wide (quasi-continuous) NegBin distributions that would consume
  too much memory with discrete RPS calculation (returning a missing
  value with a warning). [both issues spotted by F. Rousseu]

- Plots of legacy `"disProg"` and `"survRes"` objects are now generated
  via internal `disProg2sts()` conversion and `stsplot_time()`.
  This fixes their x-axis labels for the default `xaxis.years=TRUE`.
  The obsolete arguments `startyear` and `firstweek` are now ignored with
  a warning.

- The default legend of `stsplot_time1()` did not show the fill color
  in the non-default case `!is.na(col[1])`.

- Multivariate `hhh4()` with neighbourhood component treated
  `NA` counts as zero when calculating the weighted sum over units.
  A missing count at *t-1* in any unit now gives `NA` values for the
  neighbourhood terms of all units at time *t*, thus reducing `nobs()`.

## Deprecated & Defunct

- `create.disProg()` is deprecated. Methods for legacy `"disProg"` objects
  are kept for backwards compatibility, but new projects should use `sts()`.

- The long-deprecated `qlomax()` implementation has been removed.


# surveillance 1.19.1 (2021-03-30)

## Documentation

- The project website at <https://surveillance.R-Forge.R-project.org/>
  has been overhauled using [**pkgdown**](https://pkgdown.r-lib.org/).

## Bug Fixes

- The `CRS` of `data(imdepi)` and
  `data(measlesWeserEms)` have been updated via
  `sp::rebuild_CRS()` to avoid warnings when **rgdal**
  is loaded with new PROJ and GDAL libraries.

- `simEpidataCS()` now internally resets the CRS
  (temporary), which avoids spurious warnings and also reduces its
  runtime by about 25%.

- Fix encoding error in `vignette("twinstim")` for CRAN's
  non-UTF8 Linux test machine.

- This version of **surveillance** (formally) requires the new
  [**spatstat**](https://CRAN.R-project.org/package=spatstat) umbrella package to avoid collisions of old
  **spatstat** and its new sub-packages (we only use
  [**spatstat.geom**](https://CRAN.R-project.org/package=spatstat.geom)).
  The **spatstat** dependence will be dropped in the future.

- The `epoch<-` replacement method for `"sts"` objects now accepts
  a `"Date"` vector. The standard plots may give nicer x-axis annotation
  if indexed by dates. See the `xaxis.*` arguments of `stsplot_time()`.

- `tidy.sts()` (and thus `autoplot.sts()`) failed for date-indexed
  `"sts"` objects with non-standard frequencies. [spotted by Junyi Lu]


# surveillance 1.19.0 (2021-01-29)

## New Features

- The `nowcast()` function with
  `method="bayes.trunc.ddcp"` now adds support for negative
  binomial response distribution instead of Poisson. Furthermore,
  additional components of the design matrix for the discrete time
  survival model can be provided, which allows the inclusion of, e.g.,
  day of the week effects. Finally, the order of
  the polynomial created by the change-points in the discrete time
  survival model can now be specified. For further details see the
  work of Guenther et al. (2020) about nowcasting the Covid-19
  outbreak in Bavaria, Germany.

- `animate.sts()` can position the `timeplot` on
  other sides of the map.

## Minor Changes

- The weighted sum in the `ne`ighbourhood component of
  `hhh4()` models is computed more efficiently.

- `simEpidataCS()` (and thus `simulate.twinstim()`)
  uses a slightly more efficient location sampler for models with
  `siaf = siaf.constant()`.
  Simulation results will differ from previous package versions even
  if the same random `seed` is used.

- The default `main` title for `stsplot_space()` now
  uses the ISO year-week format for weekly `"sts"` data.

## Bug Fixes

- Bug fix in the `farringtonFlexible()`-function, which for
  the argument `thresholdMethod=="nbPlugin"` and
  `thresholdMethod=="muan"` unfortunately computed the limit as
  an `(1-alpha/2)` prediction interval instead of the documented
  `(1-alpha)` prediction interval. This affects four threshold
  values in Table 2 of `vignette("monitoringCounts")`.
  The default method `"delta"` worked as expected.

- In `hhh4()` models without AR component,
  the matrix of fitted values could lack column names.

- Experimental time-varying neighbourhood weights in
  `hhh4()` were indexed differently in model fitting and in the
  `simulate()` method (undocumented behaviour). Both now use the
  latter variant, where the mean at time *t* uses products of
  weights at time *t* and observed counts at time *t-1*.
  [reported by Johannes Bracher]

- For weekly `sts` indexed via `start` and `freq=52`,
  `epoch(sts, as.Date=TRUE)` now interprets the `start` week
  according to ISO 8601. For example, `start = c(2020, 5)`
  corresponds to 2020-01-27, not 2020-02-03.
  This affects `as.xts.sts()` and the time plot in
  `animate.sts()`.

- `stsplot_space()` automatically extends manual color
  breaks (`at`), if the intervals do not cover the data range.

- `simEndemicEvents()` and thus
  `epitest(..., method="simulate")` are no longer slowed down by
  intermediate `sp::CRS()` computations.

## Package Infrastructure

- Removed unused **rmapshaper** from "Suggests" and moved
  **xts** to "Enhances" (used only for `as.xts.sts`).

- Switched testing framework from (nowadays heavy)
  **testthat** to [**tinytest**](https://CRAN.R-project.org/package=tinytest). Together with moving
  **ggplot2** to "Enhances" (used only for `autoplot.sts`) ---
  and only then --- this switch further reduces the total number of
  required packages for a complete check (i.e., installing with
  `dependencies = TRUE`) in a *factory-fresh* R environment
  from 119 to 94.

- [**spatstat**](https://CRAN.R-project.org/package=spatstat) was split into several sub-packages, of which
  we only need to import [**spatstat.geom**](https://CRAN.R-project.org/package=spatstat.geom). This new package
  requires `R >= 3.5.0`, though.

- **surveillance** now requires `R >= 3.6.0`.


# surveillance 1.18.0 (2020-03-18)

## New Features

- New spatial interaction function for `twinstim()`:
  `siaf.exponential()` implements the exponential kernel
  *f(x) = exp(-x/&sigma;)*, which is a useful alternative
  if the two-parameter power-law kernel is not identifiable.

- The `plot`-type `"maps"` for `"hhh4"` fits,
  `plotHHH4_maps()`, now allows for map-specific color keys
  via `zmax = NA` (useful for `prop = TRUE`).

## Bug Fixes

- The `nowcast()`-function now also works for
  `method="bayes.trunc.ddcp"` method when the number of breakpoints
  is greater than 1.

- The `amplitudeShift` transformation for sine-cosine
  coefficient pairs in the `summary` of multivariate
  `"hhh4"` models was incorrect in the rare case that the model
  used unit-specific seasonal terms (`addSeason2formula` with
  `length(S) > 1`).

## Deprecated & Defunct

- The original `algo.hhh()` implementation of the HHH
  model has been removed from the package. The function `hhh4()`
  provides an improved and much extended implementation since 2012.


# surveillance 1.17.3 (2019-12-16)

## Bug Fixes

- The `head()`-method for `"epidataCS"` objects did
  not work with a negative `n` argument.

- Fix for `"matrix"` changes in R-devel.


# surveillance 1.17.2 (2019-11-11)

## Minor Changes

- For multivariate time series, `sts()` now checks for
  mismatches in column names of supplied matrices (`observed`,
  `population`, `neighbourhood`, ...). This is to catch
  input where the units (columns) are ordered differently in
  different slots, which would flaw subsequent analyses.

## Bug Fixes

- `simulate.twinSIR()` ignored the `atRiskY` indicator
  of the underlying `"epidata"`, so always assumed a completely
  susceptible population. Initially infectious individuals are now
  inherited. For the previous behaviour, adjust the supplied
  `data` via `data$atRiskY <- 1`.


# surveillance 1.17.1 (2019-09-13)

## New Features

- New one-parameter power-law kernel `siaf.powerlaw1()`
  with fixed `sigma = 1`. Useful if `sigma` is difficult to
  estimate with `siaf.powerlaw()`.

## Bug Fixes

- `pit()`'s default `ylab` was wrong (default are
  densities not relative frequencies).

- `R0()` for `"twinstim"` fits with specified
  `newevents` now handles levels of epidemic factor variables
  automatically via the new `xlevels` attribute stored in the
  fitted model.

- Some S3 methods for the `"sts"` class are now formally
  registered and identical to the established S4 methods.

- Minor additions and fixes in the package documentation.

## Deprecated & Defunct

- `hcl.colors()`, exported since 1.14.0, has been renamed
  `.hcl.colors()` and is now internal again, to avoid a name
  clash with R's own such function introduced in R 3.6.0.


# surveillance 1.17.0 (2019-02-22)

## New Features

- `W_powerlaw(..., from0 = TRUE)` enables more parsimonious
  `hhh4` models in that the power-law weights are modified to
  include the autoregressive (0-distance) case (see
  `vignette("hhh4_spacetime")`). The unstructured distance
  weights `W_np()` gained `from0` support as well.

- `sts()` creation can now handle `epoch` arguments of
  class `Date` directly.

- The `ranef()`-method for `"hhh4"` fits gained a
  logical argument `intercept` to extract the unit-specific
  intercepts of the log-linear predictors instead of the default
  zero-mean deviations around the fixed intercepts.
  The corresponding `plot` method (`type="ri"`) gained an
  argument `exp`: if set to `TRUE` random effects are
  `exp`-transformed and thus show multiplicative effects.
  [based on feedback by Tim Pollington]

## Minor Changes

- `W_np()`'s argument `to0` has been renamed to
  `truncate`. The old name still works but is deprecated.

- `plotHHH4_ri()` now uses `cm.colors(100)` as
  `col.regions`, and 0-centered color breaks by default.

- The help pages of `twinSIR()` and related functions now
  give examples based on `data("hagelloch")` instead of using the
  toy dataset `data("fooepidata")`. The latter is now obsolete and
  will be removed in future versions of the package.

- The elements of the `control` list stored in the result
  of `algo.farrington()` are now consistently ordered as in the
  default `control` argument.

## Bug Fixes

- Using negative indices to exclude time points from an
  `"sts"` object (e.g., `x[-1,]`) is now supported and
  equivalent to the corresponding subset expression of retained
  indexes (`x[2:nrow(x),]`) in resetting the `start` and
  `epoch` slots. [reported by Johannes Bracher]

- For weekly `"sts"` data with `epochAsDate=TRUE`,
  the `as.data.frame()` method computed `freq` by
  `"%Y"`-year instead of by `"%G"`-year, which was
  inconsistent with the `epochInPeriod` variable.

- For *non*-weekly `"sts"` data with `epochAsDate=TRUE`,
  `year()` as well as the `year` column of the
  `tidy.sts()` output corresponded to the ISO week-based year.
  It now gives the calendar year.

- `sts_creation()` hard-coded `start = c(2006, 1)`.

- `aggregate()`ing an `"sts"` object over time now
  recomputes fractions from the cumulated population values if and
  only if this is no `multinomialTS` and already contains
  population fractions. The same rule holds when subsetting units of
  an `"sts"` object. The `aggregate`-method previously
  failed to recompute fractions in some cases.

- For `farringtonFlexible()` with multivariate time series,
  only the last unit had stored the additional control items
  (exceedence scores, p-values, ...), all others were 0.
  [reported by Johannes Bracher]

- The supplementary p-values returned by `farringtonFlexible()`
  in `control$pvalue` were wrong for the default approach,
  where `thresholdMethod="delta"` (the original Farrington method)
  and a power transformation was applied to the data
  (`powertrans != "none"`). Similarly, `algo.farrington()`
  returned wrong predictive probabilities in `control$pd[,1]` if
  a power transformation was used. [reported by Lore Merdrignac]

- The `control` argument list of `algo.farrington()`
  as stated in the formal function definition was incomplete
  (`plot` was missing) and partially out of sync with the default
  values that were actually set inside the function (`b=5` and
  `alpha=0.05`). This has been fixed. Results of
  `algo.farrington()` would only be affected if the function was
  called without any `control` options (which is hardly
  possible). So this can be regarded as a documentation error.
  The formal `control` list of the `farrington()`
  wrapper function has been adjusted accordingly.

- The `control` argument lists of `farringtonFlexible()`
  and `bodaDelay()` as stated in the formal function definitions
  were partially out of sync with respect to the following default
  values that were actually set inside these functions: `b=5`
  (not 3), `alpha=0.05` (not 0.01), `pastWeeksNotIncluded=w`
  (not 26), and, for `bodaDelay()` only, `delay=FALSE` (not
  `TRUE`). This has been fixed. Results would only be affected if
  the functions were called without any `control` options (which
  is hardly possible). So this can be regarded as a documentation error.

- `pairedbinCUSUM()` did not properly subset the `sts`
  object if a `range` was specified, and forgot to store the
  `control` arguments in the result.

- `wrap.algo()` now aborts if the monitored range is
  not supplied as a numeric vector.

- In `vignette("monitoringCounts")`: several
  inconsistencies between code and output have been fixed.

- `epidataCS2sts()` no longer transfers the
  `stgrid$BLOCK` indices to the `epoch` slot of the
  resulting `"sts"` object (to avoid `epoch[1] != 1`
  scenarios).

- The `ranef()` matrix extracted from fitted `"hhh4"`
  models could have wrong column names.

## Deprecated & Defunct

- Several ancient functions deprecated in 1.16.1 are now
  defunct: `compMatrix.writeTable()`,
  `makePlot()`, `test()`, `testSim()`,
  `readData()` (the raw txt files have been removed as well),
  `correct53to52()`, `enlargeData()`, `toFileDisProg()`.


# surveillance 1.16.2 (2018-07-24)

## Minor Changes

- `autoplot.sts()` gained a `width` argument to adjust
  the bar width, which now defaults to 7 for weekly time series
  (previously was 90% of that so there were gaps between the bars).

- `"epidataCS"` generation now (again) employs
  [**spatstat**](https://CRAN.R-project.org/package=spatstat)'s `bdist.points()`, which has been
  accelerated in version 1.56-0. If you use the
  `twinstim()`-related modelling part of **surveillance**, you
  are thus advised to update your **spatstat** installation.

- The `boda()` examples in
  `vignette("monitoringCounts")` have been updated to also work
  with recent versions of **INLA**.

## Bug Fixes

- Offsets in `hhh4`'s epidemic components were ignored by
  `simulate.hhh4()` [spotted by Johannes Bracher] as well as
  in dominant eigenvalues ("maxEV").

- The color key in `fanplot()` is no longer distorted by
  `log="y"`.


# surveillance 1.16.1 (2018-05-28)

## Bug Fixes

- `autoplot.sts()` now sets the calling environment as
  the `plot_env` of the result.

- Several `twinstim`-related functions finally allow for
  prehistory events (long supported by `twinstim()` itself):
  `as.epidataCS()`, `glm_epidataCS()`,
  `as.epidata.epidataCS()`.

- The `summary()` for SI[R]S-type `"epidata"` failed
  if there were initially infectious individuals.

## Deprecated & Defunct

- Several ancient functions have been deprecated and may be
  removed in future versions of **surveillance**: `qlomax()`,
  `readData()`, `toFileDisProg()`, `correct53to52()`,
  `enlargeData()`, `compMatrix.writeTable()`,
  `test()`, `testSim()`, `makePlot()`.


# surveillance 1.16.0 (2018-01-24)

## New Features

- The `as.data.frame()` method for `"sts"` objects
  gained a `tidy` argument, which enables conversion to the
  long data format and is also available as function `tidy.sts()`.

- A [**ggplot2**](https://CRAN.R-project.org/package=ggplot2) variant of `stsplot_time()` is now
  available via `autoplot.sts()`.

- `as.epidata.data.frame()` gained an argument
  `max.time` to specify the end of the observation period (which
  by default coincides with the last observed event).

- The now exported function `fanplot()` wraps
  [**fanplot**](https://CRAN.R-project.org/package=fanplot)`::fan()`. It is used by
  `plot.oneStepAhead()` and `plot.hhh4sims()`, which now
  have an option to add the point forecasts to the fan as well.

- `plotHHH4_fitted()` (and `plotHHH4_fitted1()`)
  gained an option `total` to sum the fitted components over all
  units.

## Significant Changes

- Package [**polyCub**](https://CRAN.R-project.org/package=polyCub) is no longer automatically attached
  (only imported).

- `scores.oneStepAhead()` no longer reverses the ordering
  of the time points by default, as announced in 1.15.0.

## Minor Changes

- Some code in `vignette("monitoringCounts")` has been
  adjusted to work with the new version of [**MGLM**](https://CRAN.R-project.org/package=MGLM) (0.0.9).

- Added a `[`-method for the `"hhh4sims"` class to
  retain the attributes when subsetting simulations.

## Bug Fixes

- `aggregate(stsObj, by = "unit")` no longer results in
  empty colnames (set to `"overall"`).
  The obsolete map is dropped.

- The `subset` argument of `twinSIR()` was partially
  ignored:

  -   If `nIntervals = 1`, the model `summary()`
      reported the total number of events.

  -   Automatic `knots`, model `residuals()`, as well as
      the rug in `intensityplot()` were computed from the whole set
      of event times.

- The `as.epidata.data.frame()` converter did not actually
  allow for latent periods (via `tE.col`). This is now possible
  but considered experimental (methods for `"epidata"` currently
  ignore latent periods).

- The `all.equal()` methods for `"hhh4"` and
  `"twinstim"` objects now first check for the correct classes.


# surveillance 1.15.0 (2017-10-06)

## New Features

- `siaf.gaussian()` now also employs a `polyCub.iso()`
  integration routine by default (similar to the powerlaw-type
  kernels), instead of adaptive midpoint cubature.
  This increases precision and considerably accelerates estimation of
  `twinstim()` models with a Gaussian spatial interaction
  function. Models fitted with the new default
  (`F.adaptive=FALSE, F.method="iso"`)
  will likely differ from previous fits (`F.adaptive=TRUE`),
  and the numerical difference depends on
  the adaptive bandwidth used before (the default `adapt=0.1`
  yielded a rather rough approximation of the integral).

- Added `quantile()`, `confint()`, and `plot()`
  methods for `"oneStepAhead"` predictions.

- Exported the function `simEndemicEvents()` to simulate a
  spatio-temporal point pattern from an endemic-only
  `"twinstim"`; faster than via the general
  `simulate.twinstim()` method.

## Minor Changes

- `twinstim(..., siaf = siaf.gaussian())`
  uses a larger default initial value for the kernel's standard
  deviation (based on the size of the observation region).

- Non-default parametrizations of `siaf.gaussian()` are
  deprecated, i.e., always use `logsd=TRUE` and `density=FALSE`.

- `twinstim()` uses a smaller default initial value for the
  epidemic intercept, which usually allows for faster convergence.

- `update.hhh4()` now allows `subset.upper` values
  beyond the originally fitted time range (but still within the time
  range of the underlying `"sts"` object).

- `scores.oneStepAhead()` by default reverses the ordering
  of the time points. This awkward behaviour will change in the next
  version, so the method now warns if the default `reverse=TRUE`
  is used without explicit specification.

- Minor improvements in the documentation and some vignettes:
  corrected typos, simplified example code, documented some methods.

## Bug Fixes

- The C-routines introduced in version 1.14.0 used `==`
  comparisons on parameter values to choose among case-specific
  formulae (e.g., for *d==2* in `siaf.powerlaw()`).
  We now employ an absolute tolerance of 1e-7 (which should fix
  the failing tests on Solaris).

- Interaction functions for `twinstim()`, such as
  `siaf.powerlaw()` or `tiaf.exponential()`, no longer live in
  the global environment as this risks using masked base functions.


# surveillance 1.14.0 (2017-06-29)

## Documentation

- The replication code from Meyer et al. (2017, JSS)
  is now included as `demo("v77i11")`.
  It exemplifies the spatio-temporal endemic-epidemic modelling
  frameworks `twinstim`, `twinSIR`, and `hhh4`
  (see also the corresponding vignettes).

## New Features

- Pure C-implementations of integration routines for spatial
  interaction functions considerably accelerate the estimation of
  `twinstim()` models containing
  `siaf.powerlaw()`, `siaf.powerlawL()`, or `siaf.student()`.

- The color palette generating function used by `sts`
  plots, `hcl.colors`, is now exported.

- The utility function `clapply` (*c*onditional
  `lapply`) is now exported.

- Some utility functions for `hhh4` fits are now exported
  (`update.hhh4`, `getNEweights`, `coefW`),
  as well as several internal functions for use by `hhh4` add-on
  packages (`meanHHH`, `sizeHHH`, `decompose.hhh4`).

- The `"fan"`-type plot function for `"hhh4sims"`
  gained a `key.args` argument for an automatic color key.

- New auxiliary function `makeControl()`, which may be
  used to specify a `hhh4()` model.

## Minor Changes

- `twinstim()` now throws an informative error message when
  trying to fit a purely epidemic model to data containing endemic
  events (i.e., events without ancestors). The `help("twinstim")`
  exemplifies such a model.

## Bug Fixes

- `siaf.powerlaw()$deriv` returned `NaN` for the
  partial derivative wrt the decay parameter *d*, if *d*
  was large enough for *f* to be numerically equal to 0.
  It will now return 0 in this case.

- `twinstim()` could fail (with an error from
  `duplicated.default`) if the fitted time range was
  substantially reduced via the `T` argument.

- The `"simEpidataCSlist"` generated by
  `simulate.twinstim(..., simplify = TRUE)` was missing the
  elements `bbox` and `control.siaf`.


# surveillance 1.13.1 (2017-04-28)

## Documentation

- The paper on "Spatio-Temporal Analysis of Epidemic
  Phenomena Using the R Package **surveillance**" (by Sebastian
  Meyer, Leonhard Held, and Michael H&ouml;hle) will appear
  in the upcoming volume of the *Journal of Statistical Software*.
  The main sections 3 to 5 of the paper are contained in the package
  as `vignette("twinstim")`, `vignette("twinSIR")`, and
  `vignette("hhh4_spacetime")`, respectively.

## New Features

- The `calibrationTest()` and `pit()` methods for
  `"oneStepAhead"` forecasts gained an argument `units`
  to allow for unit-specific assessments.

- A default `scores`-method is now available to compute a
  set of proper scoring rules for Poisson or NegBin predictions.

- New plot `type = "fan"` for simulations from
  `"hhh4"` models to produce a fan chart using the
  [**fanplot**](https://CRAN.R-project.org/package=fanplot) package.

## Minor Changes

- `scores.hhh4()` sets rownames for consistency with
  `scores.oneStepAhead()`.

## Bug Fixes

- The `"Lambda.const"` matrix returned by
  `getMaxEV_season()` was wrong for models with asymmetric
  neighbourhood weights. [spotted by Johannes Bracher]\
  Dominant eigenvalues (`"maxEV"`) were not affected by this bug.


# surveillance 1.13.0 (2016-12-20)

## New Features

- `earsC` now has two new arguments thanks to Howard
  Burkom: the number of past time units to be used in calculation is
  now not always 7, it can be chosen in the `baseline` parameter.
  Furthermore, the `minSigma` parameter allows to get a threshold
  in the case of sparse data. When one doesn't give any value for those
  two parameters, the algorithm works like it used to.

- `animate.sts()` gained support for date labels in the
  bottom `timeplot`.

- `stsplot_space()` and `animate.sts()` can now
  generate incidence maps based on the population information
  stored in the supplied `"sts"` object.
  Furthermore, `animate.sts()` now supports time-varying
  population numbers.

## Minor Changes

- `hhh4()` guards against the misuse of
  `family = factor("Poisson")` for univariate time series.
  Previously, this resulted in a negative binomial model by
  definition, but is now interpreted as `family = "Poisson"`
  (with a warning).

## Bug Fixes

- `animate.sts()` now supports objects with missing values
  (with a warning). Furthermore, the automatic color breaks have been
  improved for incidence maps, also in `stsplot_space()`.

- The `as.data.frame`-method for the `"sts"` class,
  applied to classical time-index-based `"sts"` objects
  (`epochAsDate=FALSE`), ignored a `start` epoch different
  from 1 when computing the `epochInPeriod` indexes.
  Furthermore, the returned `epochInPeriod` now is a fraction of
  `freq`, for consistency with the result for objects with
  `epochAsDate=TRUE`.

- `simulate.hhh4()` did not handle shared overdispersion
  parameters correctly. The different parameters were simply recycled
  to the number of units, ignoring the factor specification from the
  model's `family`. [spotted by Johannes Bracher]

- Simulations from *endemic-only* `"hhh4"` models
  with unit-specific overdispersion parameters used wrong
  variances. [spotted by Johannes Bracher]

- `oneStepAhead()` predictions of `type`
  `"rolling"` (or `"first"`) were incorrect for time points
  `tp` (`tp[1]`) beyond the originally fitted time range
  (in that they were based on the original time range only).
  This usage of `oneStepAhead()` was never really supported and
  is now caught when checking the `tp` argument.

- `plot.hhh4simslist()` ignored its `par.settings`
  argument if `groups=NULL` (default).


# surveillance 1.12.2 (2016-11-14)

## New Features

- The internal auxiliary function, which determines the sets of
  potential source events in `"epidataCS"` has been implemented
  in C++, which accelerates `as.epidataCS()`,
  `permute.epidataCS()`, and therefore `epitest()`.
  This is only really relevant for `"epidataCS"` with a large
  number of events (>1000, say).

- Negative-binomial `hhh4()` models may not converge for
  non-overdispersed data (try, e.g.,
  `set.seed(1); hhh4(sts(rpois(104, 10)), list(family="NegBin1"))`).
  The resulting non-convergence warning message now mentions low
  overdispersion if this is detected. [suggested by Johannes Bracher]

- An additional `type="delay"` option was added to the
  `plot` method of `stsNC` objects. Furthermore, an
  `animate_nowcasts` function allows one to animate a sequence of
  nowcasts.

## Minor Changes

- In the `animate`-method for `"sts"` objects,
  the default top padding of **lattice** plots is now disabled for the
  bottom `timeplot` to reduce the space between the panels.
  Furthermore, the new option `fill` can be used to make the
  panel of the `timeplot` as large as possible.

## Bug Fixes

- `bodaDelay()`: fixed spurious warnings from `rnbinom()`.

- `vignette("monitoringCounts")`: fixed `boda`-related
  code and cache to obtain same results as in corresponding JSS paper.


# surveillance 1.12.1 (2016-05-18)

## Documentation

- The new `vignette("monitoringCounts")` illustrates the
  monitoring of count time series in R with a particular focus on
  aberration detection in public health surveillance.
  This vignette corresponds to a recently accepted manuscript
  for the *Journal of Statistical Software*
  (Salmon, Schumacher, and H&ouml;hle, 2016).

## Minor Changes

- Non-convergent `hhh4()` fits now obey the structure of
  standard `"hhh4"` objects. In particular, such fits now also
  contain the `control` and `stsObj` elements, allowing
  for model `update()`s of non-convergent fits.

- `knox()` warns about symmetric input matrices.

## Bug Fixes

- The code of `boda()` (with `samplingMethod="joint"`)
  and `bodaDelay()` (with `inferenceMethod="INLA"`)
  has been adjusted to a change of arguments of **INLA**'s
  `inla.posterior.sample` function. Accordingly, the minimum
  **INLA** version required to run `boda()` and
  `bodaDelay()` is 0.0-1458166556.

- The functions returned by `W_powerlaw()` now have the
  package namespace as their environment to support situations where
  the package is not attached.

- Attaching package [**nlme**](https://CRAN.R-project.org/package=nlme) after **surveillance** no
  longer masks `"hhh4"`'s `ranef`-method. (We now import the
  `fixef` and `ranef` generics from **nlme**.)


# surveillance 1.12.0 (2016-04-02)

## Documentation

- Several new vignettes illustrate *endemic-epidemic*
  modeling frameworks for spatio-temporal surveillance data:

  `vignette("twinstim")`

  :   describes a spatio-temporal
      point process regression model.

  `vignette("twinSIR")`

  :   describes a multivariate
      temporal point process regression model.

  `vignette("hhh4_spacetime")`

  :   describes an areal
      time-series model for infectious disease counts.

  These vignettes are based on a recently accepted manuscript
  for the *Journal of Statistical Software*
  (Meyer, Held, and H&ouml;hle, 2016).

- Improved the documentation on various help pages.

- The `hhh4()`-based analysis of `data("fluBYBW")`
  has been moved to a separate demo script 'fluBYBW.R'. Due to
  the abundance of models and the relatively long runtime, we
  recommend to open the script in an editor rather than running
  all the code at once using `demo("fluBYBW")`.

## New Features

- Overhaul of the `"sts"` implementation. This mostly
  affects package-internal code, which is simpler, cleaner and better
  tested now, but requires R >= 3.2.0 (due to `callNextMethod()`
  bugs in older versions of R).
  Beyond that, the user-level constructor function `sts()` now
  has explicit arguments for clarity and convenience.
  For instance, its first argument sets the `observed` slot and
  no longer needs to be named, i.e.,
  `sts(mycounts, start=c(2016,3), frequency=12)`
  works just like for the classical `ts()` function.

- `stsplot_time(..., as.one=TRUE)` is now implemented
  (yielding a simple `matplot` of multiple time series).

## Minor Changes

- `plotHHH4_season()` now by default draws a horizontal
  reference line at unity if the multiplicative effect of component
  seasonality is shown (i.e., if `intercept=FALSE`).

- Since **surveillance** 1.8-0, `hhh4()` results are of
  class `"hhh4"` instead of `"ah4"` (renamed).
  Legacy methods for the old class name `"ah4"` have been removed.

- The internal model preparation in `twinstim()` is more
  efficient (the distance matrix of the events is only computed if
  event sources actually need to be updated).

## Bug Fixes

- `stsplot_spacetime()` now recognizes its `opts.col`
  argument.

- Conversion from `"ts"` to `"sts"` using
  `as(ts, "sts")` could set a wrong start time. For instance,
  `as(ts(1:10, start=c(1959,2), frequency=4), "sts")@start` was
  `c(1959,1)`.

- `algo.twins()` now also accepts `"sts"` input and
  the automatic legend in the first plot of `plot.atwins()` works
  again.

- The experimental `profile`-method for `"twinstim"`
  objects did not work if embedded `twinstim()` fits issued warnings.


# surveillance 1.11.0 (2016-02-08)

## New Features

- `update.epidata()` can now handle a distance matrix
  `D` in the form of a classed `"Matrix"`.
  [suggested by George Wood]

- `glrnb()` can now handle `ret="cases"` for the
  generalized likelihood ratio detector based on the negative binomial
  distribution. It's based on a brute-force search and hence might be
  slow in some situations.

- `boda()` and `bodaDelay()` now support an
  alternative method (`quantileMethod="MM"`) to compute quantiles
  based on the posterior distribution. The new method samples
  parameters from the posterior distribution and then computes the
  quantile of the mixture distribution using bisectionning, which is
  faster and yields similar results compared to the original method
  (`quantileMethod="MC"`, still the default).

## Minor Changes

- Revised `vignette("hhh4")`, updated the package
  description as well as some references in the documentation.
  Also updated (the cache of) the slightly outdated
  `vignette("surveillance")` to account for the corrected version
  of `algo.bayes()` implemented since **surveillance** 1.10-0.

## Bug Fixes

- Fixed bug in `categoricalCUSUM()`, which ignored alarms
  generated for the last time point in `range`. Furthermore, the
  exact computation in case of returns of the type `"value"` for
  the binomial are now checked through an attribute.

- Fixed bug in the `estimateGLRNbHook` function of
  `algo.glrnb`, which ignored potential fixed `alpha`
  values. If `alpha` is fixed this is now taken into consideration
  while fitting the negative binomial function. See revised help files
  for the details.

- Made a hot-fix such that the `algo.quality` function now
  also works for `sts` objects and if the `state` or
  `alarm` slots consists of TRUE/FALSE instead of 0/1.

- `intensity.twinstim()` did not work for non-endemic models.

- A parallelized `epitest()` could fail with a strange
  error message if some replications were left unassigned.
  This seems to happen if forking is used (`mclapply`) with
  insufficient memory. Incomplete replications are now ignored with a
  warning.


# surveillance 1.10-0 (2015-11-04)

This package is now maintained by Sebastian Meyer, who has been an active
co-author since version 1.3. We thank Michael H&ouml;hle for 10 years of
maintenance ever since he created **surveillance** and published the
package on CRAN in November 2005 for R 2.2.0.

## New Features

- Calibration tests for count data (Wei and Held, 2014, Test)
  are now implemented and available as `calibrationTest()`.
  In addition to a default method taking pure counts and predictive
  means and dispersion parameters, there are convenient methods for
  `"hhh4"` and `"oneStepAhead"` objects.

- Shared overdispersion across units in negative binomial
  `hhh4()` time series models (by specifying a factor variable
  as the `family` argument).

- `scores()` and `pit()` are now generic and have convenient
  methods for `"oneStepAhead"` predictions and `"hhh4"` fits.

- The initial values used for model updates during the
  `oneStepAhead()` procedure can now be specified directly
  through the `which.start` argument (as an alternative to the
  previous options `"current"` and `"final"`).

- `plotHHH4_fitted()` (and `plotHHH4_fitted1()`)
  gained an option `decompose` to plot the contributions from
  each single unit (and the endemic part) instead of the default
  endemic + AR + neighbours decomposition.
  Furthermore, a formatted time axis similar to `stsplot_time1()`
  can now be enabled via the new argument `xaxis`.

- The new `plot` `type` `"maps"` for `"hhh4"`
  fits shows maps of the fitted mean components averaged over time.

- New `plot`-method for simulations from `"hhh4"`
  models (using `simulate.hhh4(..., simplify = TRUE)`, which now
  has a dedicated class: `"hhh4sims"`) to show the final
  size distribution or the simulated time series
  (possibly stratified by groups of units).
  There is also a new `scores`-method to compute proper scoring
  rules based on such simulations.

- The argument `idx2Exp` of `coef.hhh4()` may now be
  conveniently set to `TRUE` to exp-transform all coefficients.

- Added a `coeflist()`-method for `"hhh4"` fits.

- The generator function `sts()` can now be used to
  initialize objects of class `"sts"` (instead of writing
  `new("sts", ...)`).

- Additional arguments of `layout.scalebar()` now allow to
  change the style of the labels.

- A pre-computed distance matrix `D` can now be used as
  input for the `as.epidata()` converter -- offering an alternative
  to the default Euclidean distance based on the individuals coordinates.
  (Request of George Wood to support `twinSIR` models on networks.)

## Minor Changes

- The first argument of `scores()` is now called `x`
  instead of `object` (for consistency with `calibrationTest()`).

- The result of `oneStepAhead()` now has the dedicated
  class attribute `"oneStepAhead"` (previously was just a list).

- Changed interpretation of the `col` argument of
  `plotHHH4_fitted()` and `plotHHH4_fitted1()` (moved color
  of "observed" to separate argument `pt.col` and reversed
  remaining colors). The old `col` specification as a vector of
  length 4 still works (caught internally) but is undocumented.

- The `epoch` slot of class `"sts"` is now initialized to
  `1:nrow(observed)` by default and thus no longer needs to be
  explicitly set when creating a `new("sts", ...)` for this
  standard case.

- Initialization of `new("sts", ...)` now supports the
  argument `frequency` (for consistency with `ts()`).
  Note that `freq` still works (via partial argument matching)
  and that the corresponding `"sts"` slot is still called `freq`.

- If `missing(legend.opts)` in `stsplot_time1()`, the
  default legend will only be produced if the `"sts"` object
  contains information on outbreaks, alarms, or upperbounds.

- The default `summary()` of a `"twinstim"` fit is more
  concise since it no longer includes the number of log-likelihood and
  score function evaluations and the elapsed time during model fitting.
  Set the new `runtime` argument of `summary.twinstim()` to
  `TRUE` to add this information to the summary as before.

- The `animate`-method for `"sts"` objects gained an
  argument `draw` (to disable the default instantaneous plotting)
  and now invisibly returns the sequential plot objects (of class
  `"gtable"` or `"trellis"`) in a list for post-processing.

- The flexible time axis configurations for `"sts"` plots
  introduced in version 1.8-0 now also work for classical `"sts"`
  objects with integer epochs and standard frequencies
  (try `plot(..., epochsAsDate = TRUE)`).

- `stsplot_time()` initiates `par` settings only
  if the `par.list` argument is a list.

- The new `all.equal()` method for class `"hhh4"`
  compares two fits ignoring their `"runtime"` and `"call"`
  elements (at least).

## Bug Fixes

- Fixed a bug in `algo.bayes`, where an alarm was already
  sounded if the current observation was equal to the quantile of the
  predictive posterior. This was changed in order to get *alarm_t
  = I(obs_t > quantile_t)* which is consistent with the use in
  `boda` and `bodaDelay`.

- Fixed bug in `algo.outbreakP` causing a halt in the
  computations of `value="cases"` when
  `calc.outbreakP.statistic` returned `NaN`. Now, a
  `NaN` is returned.

- `wrap.algo` argument `control.hook` used
  `control` argument defined outside it's scope (and not the one
  provided to the function). It is now added as additional 2nd argument
  to the `control.hook` function.

- `stsplot_time()` did not account for the optional
  `units` argument for multivariate `"sts"` objects
  when choosing a suitable value for `par("mfrow")`.

- `hhh4()` could have used a function `dpois()` or
  `dnbinom()` from the global environment instead of the
  respective function from package **stats**.

- The default time variable `t` created as part of the
  `data` argument in `hhh4()` was incompatible with
  `"sts"` objects having `epochAsDate=TRUE`.

- A consistency check in `as.epidata.default()` failed for
  SI-type data (and, more generally, for all data which ended with an
  I-event in the last time block). [spotted by George Wood]


# surveillance 1.9-1 (2015-06-12)

- This is a quick patch release to make the test suite run
  smoothly on CRAN's Windows and Solaris Sparc systems.

- The new `hhh4()` option to scale neighbourhood weights
  did not work for parametric weights with more than one parameter
  if `normalize=FALSE`.


# surveillance 1.9-0 (2015-06-09)

## New Features

- New functions and data for Bayesian outbreak detection in the
  presence of reporting delays (Salmon et al., 2015):
  `bodaDelay()`, `sts_observation()`, and `sts_creation()`.

- New functions implementing tests for space-time interaction:

  -   `knox()` supports both the Poisson approximation and a
      Monte Carlo permutation approach to determine the p-value,

  -   `stKtest()` wraps space-time K-function methods from
      package [**splancs**](https://CRAN.R-project.org/package=splancs) for use with `"epidataCS"`,

  -   and `epitest()` for `twinstim` models
      (makes use of the new auxiliary function `simpleR0()`).

- New function `plapply()`: a parallel and verbose version
  of `lapply()` wrapping around both `mclapply()` and
  `parLapply()` of package **parallel**.

- New converter `as.xts.sts()` to transform `"sts"`
  objects to the quasi standard `"xts"` class, e.g., to make use
  of package [**dygraphs**](https://CRAN.R-project.org/package=dygraphs) for interactive time series plots.

- New options for scaling and normalization of neighbourhood
  weights in `hhh4()` models.

- New auxiliary function `layout.scalebar()` for use as part
  of `sp.layout` in `sp::spplot()` or in the traditional
  graphics system.

### New features for `"epidataCS"`

- New argument `by` for `plot.epidataCS()`, which
  defines a stratifying variable for the events (default is the event
  type as before). It can also be set to `NULL` to make the plot
  not distinguish between event types.

- The spatial plot of `"epidataCS"` gained the arguments
  `tiles`, `pop` and `sp.layout`, and can now produce
  an `sp::spplot()` with the tile-specific population levels behind
  the point pattern.

- New function `permute.epidataCS()` to randomly permute
  time points or locations of the events (holding other marks fixed).

### New features for `twinstim()`

- New S3-generic `coeflist()` to list model coefficients by
  component. It currently has a default method and one for
  `"twinstim"` and `"simEpidataCS"`.

- New argument `newcoef` for `simulate.twinstim()` to
  customize the model parameters used for the simulation.

- New argument `epilink` for `twinstim()`, offering
  experimental support for an identity link for the epidemic
  predictor. The default remains `epilink = "log"`.

- Simulation from `"twinstim"` models and generation of
  `"epidataCS"` is slightly faster now (faster **spatstat**
  functions are used to determine the distance of events to the border).

- New option `scaled = "standardized"` in `iafplot()`
  to plot *f(x) / f(0)* or *g(t) / g(0)*, respectively.

## Minor Changes

- Initial data processing in `twinstim()` is faster
  since event sources are only re-determined if there is effective
  need for an update (due to subsetting or a change of
  `qmatrix`).

- `formatPval()` disables `scientific` notation by default.

- The `"time"` plot for `"epidataCS"` uses the
  temporal grid points as the default histogram `breaks`.

- The special `fe()` function which sets up fixed effects
  in `hhh4()` models gained an argument `unitSpecific` as a
  convenient shortcut for `which = rep(TRUE, nUnits)`.

- The convenient `plot` option of `permutationTest()`
  uses [**MASS**](https://CRAN.R-project.org/package=MASS)::`truehist()` instead of `hist()` and
  accepts graphical parameters to customize the histogram.

## Bug Fixes

- The `bodaFit` function did not draw samples from the
  joint posterior. Instead draws were from the respective posterior
  marginals. A new argument `samplingMethod` is now introduced
  defaulting to the proper 'joint'. For backwards compatibility use
  the value 'marginal'.

- The functions `as.epidataCS()` and `simEpidataCS()`
  could throw inappropriate warnings when checking polygon areas
  (only if `W` or `tiles`, respectively, contained holes).

- Non-convergent endemic-only `twinstim` models
  produced an error. [spotted by Bing Zhang]

- The `"owin"`-method of `intersectPolyCircle` could
  have returned a rectangle-type `"owin"` instead of a polygon.

- An error occurred in `twinstim()` if `finetune=TRUE`
  or choosing `optim()` instead of the default `nlminb()`
  optimizer without supplying a `control` list in `optim.args`.

- The `"time"` plot for `"epidataCS"` did not
  necessarily use the same histogram `breaks` for all strata.

- Specifying a step function of interaction via a numeric vector
  of knots did not work in `twinstim()`.

- `plot.hhh4()` did not support an unnamed `type`
  argument such as `plot(x, "season")`.

- `simEpidataCS()` did not work if `t0` was in the
  last block of `stgrid` (thus it did not work for single-cell
  grids), and mislabeled the `start` column copied to
  `events` if there were no covariates in `stgrid`.

- Evaluating `intensity.twinstim()$hFUN()` at time points
  before `t0` was an error. The function now returns
  `NA_real_` as for time points beyond `T`.

- Truncated, normalized power-law weights for `hhh4()`
  models, i.e., `W_powerlaw(maxlag = M, normalize = TRUE)`
  with `M < max(neighbourhood(stsObj))`, had wrong derivatives
  and thus failed to converge.

- `update.hhh4(..., use.estimates = TRUE)` did not
  use the estimated weight function parameters as initial values for
  the new fit. It does so now iff the weight function
  `ne$weights` is left unchanged.


# surveillance 1.8-3 (2015-01-05)

- Accommodate a new note given by R-devel checks, and set the new
  INLA additional repository in the 'DESCRIPTION' file.

- Made `linelist2sts()` work for quarters by adding extra
  `"%q"` formatting in `formatDate()`.


# surveillance 1.8-2 (2014-12-16)

## Minor Changes for `hhh4()`

- In the coefficient vector resulting from a `hhh4` fit,
  random intercepts are now named.

- Parameter `start` values in `hhh4()` are now
  matched by name but need not be complete in that case (default
  initial values are used for unspecified parameters).

- The `update.hhh4()`-method now by default does
  `use.estimates` from the previous fit. This reduces the
  number of iterations during model fitting but may lead to slightly
  different parameter estimates (within a tolerance of `1e-5`).
  Setting `use.estimates = FALSE` means to re-use the previous
  start specification.

## Minor Changes for the `"sts"` Class

- For univariate `"sts"` objects, the (meaningless)
  "head of neighbourhood" is no longer `show`n.

- The `"sts"` class now has a `dimnames`-method
  instead of a `colnames`-method. Furthermore, the redundant
  `nrow` and `ncol` methods have been removed
  (the `dim`-method is sufficient).

- If a `map` is provided when `initialize()`ing an
  `"sts"` object, it is now verified that all `observed`
  regions are part of the `map` (matched by `row.names`).

- In `stsplot_space()`, extra (unobserved) regions of the
  `map` are no longer dropped but shown with a dashed border by
  default.


# surveillance 1.8-1 (2014-10-29)

## New Features

- The `R0`-method for `"twinstim"` gained an argument
  `newcoef` to simplify computation of reproduction numbers with
  a different parameter vector (also used for Monte Carlo CI's).

- New plot `type="neweights"` for `"hhh4"` fits.

- The `scores()` function allows the selection of multiple
  `units` (by index or name) for which to compute (averaged) proper
  scores. Furthermore, one can now select `which` scores to compute.

- Added a `formula`-method for `"hhh4"` fits to
  extract the `f` specifications of the three components from the
  control list.

- The `update()`-method for fitted `"hhh4"` models
  gained an argument `S` for convenient modification of component
  seasonality using `addSeason2formula()`.

- The new auxiliary function `layout.labels()` generates an
  `sp.layout` item for `sp::spplot()` in order to draw labels.

- When generating the `pit()` histogram with a single
  predictive CDF `pdistr`, the `...` arguments can now be
  `x`-specific and are recycled if necessary using `mapply()`.
  If `pdistr` is a list of CDFs, `pit()` no longer requires
  the functions to be vectorized.

- New method `as.epidata.data.frame()`, which constructs the
  start/stop SIR event history format from a simple individual-based
  data frame (e.g., `hagelloch.df`).

- New argument `w` in `as.epidata.default()` to
  generate covariate-based weights for the force of infection in
  `twinSIR`. The `f` argument is for distance-based weights.

- The result of `profile.twinSIR()` gained a class and an
  associated `plot`-method.

## Significant Changes

- For multivariate `oneStepAhead()` predictions,
  `scores(..., individual=TRUE)` now returns a 3d array instead
  of a collapsed matrix. Furthermore, the scores computed by default
  are `c("logs","rps","dss","ses")`, excluding the normalized
  squared error score `"nses"` which is improper.

- The plot-`type="season"` for `"hhh4"` fits now by
  default plots the multiplicative effect of seasonality on the
  respective component (new argument `intercept=FALSE`).
  The default set of components to plot has also changed.

- When `as.epidata()` and `simEpidata()` calculate
  distance-based epidemic weights from the `f` functions, they no
  longer set the distance of an infectious individual to itself
  artificially to `Inf`.
  This changes the corresponding columns in the `"epidata"` in
  rows of currently infectious individuals, but the `twinSIR`
  model itself is invariant, since only rows with `atRiskY=1`
  contribute to the likelihood.

- Several modifications and corrections in `data("hagelloch")`.

## Minor Changes

- Better plotting of `stsNC` objects by writing an own plot
  method for them. Prediction intervals are now shown jointly with the
  point estimate.

- Reduced package size by applying `tools::resaveRdaFiles`
  to some large datasets and by building the package with
  `--compact-vignettes=both`, i.e., using additional GhostScript
  compression with ebook quality, see `?tools::compactPDF`.

- Added `units` argument to `stsplot_time` to select
  only a subset of the multivariate time series for plotting.

- The `untie`-method for class `"epidataCS"` gained an
  argument `verbose` which is now `FALSE` by default.

- `"epidataCS"` objects store the `clipper` used
  during generation as attribute of `$events$.influenceRegion`.

- In `plotHHH4_fitted()`, the argument `legend.observed`
  now defaults to `FALSE`.

- The default weights for the spatio-temporal component in
  `hhh4` models now are `neighbourhood(stsObj) == 1`.
  The previous default `neighbourhood(stsObj)` does not make
  sense for the newly supported `nbOrder` neighbourhood matrices
  (shortest-path distances). The new default makes no difference for
  (old) models with binary adjacency matrices in the neighbourhood
  slot of the `stsObj`.

- The default for nonparametric weights `W_np()` in
  `hhh4()` is now to assume zero weight for neighbourhood orders
  above `maxlag`, i.e., `W_np()`'s argument `to0` now
  defaults to `TRUE`.

- Added a `verbose` argument to `permutationTest()`,
  which defaults to `FALSE`. The previous behaviour corresponds
  to `verbose=TRUE`.

- `simulate.twinstim()` now by default uses the original
  `data$W` as observation region.

- The `data("measlesWeserEms")` contain two additional
  variables in the `@map@data` slot: `"vaccdoc.2004"` and
  `"vacc1.2004"`.

- The plot-method for `"epidata"` objects now uses colored
  lines by default.

- The **surveillance** package now depends on R >= 3.0.2,
  which, effectively, is the minimum version required since
  **surveillance** 1.7-0.

- The two diagnostic plots of `checkResidualProcess()` are
  now by default plotted side by side (`mfrow=c(1,2)`) instead of
  one below the other.

## Bug Fixes

- In `farringtonFlexible` alarms are now for
  `observed>upperbound` and not for `observed>=upperbound`
  which was not correct.

- Fixed duplicate `"functions"` element resulting from
  `update.twinstim(*,model=TRUE)` and ensured that
  `"twinstim"` objects always have the same components (some may
  be `NULL`).

- `animate.epidata` works again with the
  [**animation**](https://CRAN.R-project.org/package=animation) package (`ani.options("outdir")` was
  removed in version 2.3)

- For `hhh4` models with random effects, `confint()`
  only worked if argument `parm` was specified.

- Computing one-sided AIC weights by simulation for `twinSIR`
  models with more than 2 epidemic covariates now is more robust (by
  rescaling the objective function for the quadratic programming
  solver) and twice as fast (due to code optimization).

- `simulate.twinstim(..., rmarks=NULL)` can now handle the
  case where `data` has no events within the simulation period
  (by sampling marks from all of `data$events`).

- The `lambda.h` values of simulated events in
  `"simEpidataCS"` objects were wrong if the model contained an
  endemic intercept (which is usually the case).

- Automatic choice of color breaks in the `animate`-method
  for class `"sts"` now also works for incidence maps (i.e., with
  a `population` argument).

- `hhh4()` did not allow the use of nonparametric
  neighbourhood weights `W_np()` with `maxlag=2`.

- `scores()` did not work for multivariate `oneStepAhead()`
  predictions if both `individual=TRUE` and `sign=TRUE`, and
  it could not handle a `oneStepAhead()` prediction of only one
  time point. Furthermore, the `"sign"` column of
  `scores(..., sign=TRUE)` was wrong (reversed).

- For `"epidataCS"` with only one event,
  `epidataCSplot_space()` did not draw the point.

- The trivial (identity) call
  `aggregate(stsObj, nfreq=stsObj@freq)` did not work.


# surveillance 1.8-0 (2014-06-16)

## Package Infrastructure

- Package **surveillance** now depends on newer versions of
  packages [**sp**](https://CRAN.R-project.org/package=sp) (>= 1.0-15), [**polyCub**](https://CRAN.R-project.org/package=polyCub) (>= 0.4-2),
  and [**spatstat**](https://CRAN.R-project.org/package=spatstat) (>= 1.36-0).
  The R packages **INLA** and [**runjags**](https://CRAN.R-project.org/package=runjags) are now suggested
  to support a new outbreak detection algorithm (`boda()`) and
  the new `nowcast()`ing procedure, respectively.
  The R packages for [**lattice**](https://CRAN.R-project.org/package=lattice), [**grid**](https://CRAN.R-project.org/package=grid),
  [**gridExtra**](https://CRAN.R-project.org/package=gridExtra), and [**scales**](https://CRAN.R-project.org/package=scales) are suggested for
  added visualization facilities.

- More tests have been implemented to ensure package integrity.
  We now use [**testthat**](https://CRAN.R-project.org/package=testthat) instead of the outdated package
  [**RUnit**](https://CRAN.R-project.org/package=RUnit).

- `hhh4()` fits now have class `"hhh4"` instead of
  `"ah4"`, for consistency with `twinstim()`,
  `twinSIR()`, and to follow the common convention (cp. `lm()`).
  Standard S3-methods for the old `"ah4"` name are still
  available for backwards compatibility but may be removed in the
  future.

- Plot variants for `"sts"` objects have been cleaned up:
  The functions implementing the various plot types
  (`stsplot_*`, previously named `plot.sts.*`)
  are now exported and documented separately.

## New Features

- The `nowcast` procedure has been completely re-written to
  handle the inherit right-truncation of reporting data (best
  visualized as a reporting triangle). The new code implements the
  generalized-Dirichlet and the hierarchical Bayesian approach described in
  H&ouml;hle and an der Heiden (2014). No backwards
  compatibility to the old nowcasting procedure is given.

- The package contains a new monitoring function
  `boda`. This is a first experimental surveillance
  implementation of the Bayesian Outbreak Detection Algorithm (BODA)
  proposed in Manitz and H&ouml;hle (2012). The function
  relies on the non-CRAN package **INLA**, which has to be installed first
  in order to use this function. Expect initial problems.

- New `toLatex`-method for `"sts"` objects.

- The new function `stsplot_space()` provides an improved
  map plot of disease incidence for `"sts"` objects aggregated
  over time. It corresponds to the new `type = observed ~ unit`
  of the `stsplot`-method, and supersedes
  `type = observed ~ 1|unit` (except for alarm shading).

- An `animate()`-method for the `"sts"` class provides
  a new implementation for animated maps (superseding the `plot`
  `type=observed ~ 1 | unit * time`) with an optional evolving
  time series plot below the map.

- The `plot()` method for `"sts"` objects with epochs as
  dates is now made more flexible by introducing the arguments
  `xaxis.tickFreq`, `xaxis.labelFreq` and
  `xaxis.labelFormat`. These allow the specification of
  tick-marks and labelling based on `strftime` compatible
  conversion codes -- independently if data are daily, weekly, monthly,
  etc. As a consequence, the old argument `xaxis.years` is
  removed. See `stsplot_time()` for more information.

- Inference for neighbourhood weights in `hhh4()` models:
  `W_powerlaw()` and `W_np()` both implement weights
  depending on the order of neighbourhood between regions, a power-law
  decay and nonparametric weights, i.e., unconstrained estimation of
  individual weights for each neighbourhood order.

- `hhh4()` now allows the inclusion of multiplicative
  offsets also in the epidemic components `"ar"` and `"ne"`.

- `hhh4()` now has support for `lag != 1` in the
  autoregressive and neighbor-driven components. The applied lags are
  stored as component `"lags"` of the return value (previously
  there was an unused component `"lag"` which was always 1 and
  has been removed now).

- `oneStepAhead()`:

  -   Added support for parallel computation of predictions using
      `parallel::mclapply()`.

  -   New argument `type` with a new `type`
      `"first"` to base all subsequent one-step-ahead predictions
      on a single initial fit.

  -   Nicer interpretation of `verbose` levels, and
      `txtProgressBar()`.

- The `plot()`-method for fitted `hhh4()` objects now
  offers three additional types of plots: component seasonality,
  seasonal or time course of the dominant eigenvalue, and maps
  of estimated random intercepts. It is documented and more customizable.
  Note that argument order and some names have changed:
  `i` -> `units`, `title` -> `names`.

- (Deviance) `residuals()`-method for fitted `hhh4()`
  models.

- Added methods of `vcov()` and `nobs()`
  for the `"hhh4"` class. For `AIC()` and `BIC()`, the
  default methods work smoothly now (due to changes to
  `logLik.hhh4()` documented below).

- New predefined interaction functions for `twinstim()`:
  `siaf.student()` implements a *t*-kernel for the distance
  decay, and `siaf.step()` and `tiaf.step()` provide step
  function kernels (which may also be invoked by specifying the
  vector of knots as the `siaf` or `tiaf` argument in
  `twinstim`).

- Numerical integration over polygonal domains in the `F`
  and `Deriv` components of `siaf.powerlaw()` and
  `siaf.powerlawL()` is much faster and more accurate now since
  we use the new `polyCub.iso()` instead of `polyCub.SV()`
  from package [**polyCub**](https://CRAN.R-project.org/package=polyCub).

- New `as.stepfun()`-method for `"epidataCS"` objects.

- `plot.epidataCS()`:

  -   The spatial plot has new arguments to automatically add
      legends to the plot: `legend.types` and `legend.counts`.
      It also gained an `add` argument.

  -   The temporal plot now supports type-specific sub-histograms,
      additional lines for the cumulative number of events, and an
      automatic legend.

- The new function `glm_epidataCS()` can be used to fit
  an endemic-only `twinstim()` via `glm()`.
  This is mainly provided for testing purposes since wrapping into
  `glm` usually takes longer.

## Significant Changes

- Fitted `hhh4()` objects no longer contain the associated
  `"sts"` data twice: it is now only stored as `$stsObj`
  component, the hidden duplicate in `$control$data$.sts` was
  dropped, which makes fitted objects substantially smaller.

- `logLik.hhh4()` always returns an object of class
  `"logLik"` now; for random effects models, its `"df"`
  attribute is `NA_real_`. Furthermore, for non-convergent fits,
  `logLik.hhh4()` gives a warning and returns `NA_real_`;
  previously, an error was thrown in this case.

- `oneStepAhead()`:

  -   Default of `tp[2]` is now the penultimate time point of
      the fitted subset (not of the whole `stsObj`).

  -   `+1` on rownames of `$pred` (now the same as for
      `$observed`).

- The optional `"twinstim"` result components
  `fisherinfo`, `tau`, and `functions` are always
  included. They are set to `NULL` if they are not applicable
  instead of missing completely (as before), such that all
  `"twinstim"` objects have the same list structure.

- `iafplot()` ...

  -   invisibly returns a matrix containing the plotted
      values of the (scaled) interaction function (and the confidence
      interval as an attribute).
      Previously, nothing (`NULL`) was returned.

  -   detects a type-specific interaction function and by default
      uses `types=1` if it is not type-specific.

  -   has better default axis ranges.

  -   adapts to the new step function kernels (with new arguments
      `verticals` and `do.points`).

  -   supports logarithmic axes (via new `log` argument
      passed on to `plot.default`).

  -   optionally respects `eps.s` and `eps.t`,
      respectively (by the new argument `truncated`).

  -   now uses `scaled=TRUE` by default.

- The argument `colTypes` of
  `plot.epidataCS(,aggregate="space")` is deprecated (use
  `points.args$col` instead).

- The events in an `"epidataCS"` object no longer have
  a reserved `"ID"` column.

## Minor Changes

- `hhh4()` now stores the runtime just like `twinstim()`.

- Take `verbose=FALSE` in `hhh4()` more seriously.

- `hhh4()` issues a `warning()` if non-convergent.

- The following components of a `hhh4()` fit now have names:
  `"se"`, `"cov"`, `"Sigma"`.

- The new default for `pit()` is to produce the plot.

- The `twinstim()` argument `cumCIF` now defaults to
  `FALSE`.

- `update.twinstim()` no longer uses recursive
  `modifyList()` for the `control.siaf` argument. Instead,
  the supplied new list elements (`"F"`, `"Deriv"`)
  completely replace the respective elements from the original
  `control.siaf` specification.

- `siaf.lomax()` is now defunct (it has been deprecated
  since version 1.5-2); use `siaf.powerlaw()` instead.

- Allow the default `adapt`ive bandwidth to be specified via the
  `F.adaptive` argument in `siaf.gaussian()`.

- Unsupported options (`logpars=FALSE`,
  `effRangeProb`) have been dropped from `siaf.powerlaw()`
  and `siaf.powerlawL()`.

- More rigorous checking of `tiles` in
  `simulate.twinstim()` and `intensityplot.twinstim()`.

- `as.epidataCS()` gained a `verbose` argument.

- `animate.epidataCS()` now by default does not draw
  influence regions (`col.influence=NULL`), is `verbose` if
  `interactive()`, and ignores `sleep` on non-interactive
  devices.

- The `multiplicity`-generic and its default method have
  been integrated into [**spatstat**](https://CRAN.R-project.org/package=spatstat) and are imported from there.

## Data

- The polygon representation of Germany's districts (
  `system.file("shapes", "districtsD.RData", package="surveillance")`
  ) has been simplified further. The union of `districtsD` is
  used as observation window `W` in `data("imdepi")`. The
  exemplary `twinstim()` fit `data("imdepifit")` has been
  updated as well. Furthermore, `row.names(imdepi$events)` have
  been reset (chronological index), and numerical differences
  in `imdepi$events$.influenceRegion` are due to changes in
  [**polyclip**](https://CRAN.R-project.org/package=polyclip) 1.3-0.

- The Campylobacteriosis data set `campyDE`, where absolute
  humidity is used as concurrent covariate to adjust the outbreak
  detection is added to the package to exemplify `boda()`.

- New `data("measlesWeserEms")` (of class `"sts"`),
  a corrected version of `data("measles.weser")` (of the old
  `"disProg"` class).

## Bug Fixes

- Fixed a bug in `LRCUSUM.runlength` where computations
  were erroneously always done under the in-control parameter
  `mu0` instead of `mu`.

- Fixed a bug during alarm plots (`stsplot_alarm()`),
  where the use of `alarm.symbol` was ignored.

- Fixed a bug in `algo.glrnb` where the overdispersion
  parameter `alpha` from the automatically fitted `glm.nb`
  model (fitted by `estimateGLRNbHook`) was incorrectly taken as
  `mod[[1]]$theta` instead of `1/mod[[1]]$theta`. The error is
  due to a different parametrization of the negative binomial
  distribution compared to the parametrization in H&ouml;hle
  and Paul (2008).

- The score function of `hhh4()` was wrong when fitting
  endemic-only models to a `subset` including the first time
  point. This led to "false convergence".

- `twinstim()` did not work without an endemic offset if
  `is.null(optim.args$par)`.


# surveillance 1.7-0 (2013-11-19)

## Package Infrastructure

- Package [**gpclib**](https://CRAN.R-project.org/package=gpclib) is no longer necessary for the
  construction of `"epidataCS"`-objects. Instead, we make use of
  the new dedicated package [**polyclip**](https://CRAN.R-project.org/package=polyclip) (licensed under the
  BSL) for polygon clipping operations (via
  `spatstat::intersect.owin()`). This results in a slightly
  different `$events$.influenceRegion` component of
  `"epidataCS"` objects, one reason being that
  **polyclip** uses integer arithmetic.
  Change of `twinstim()` estimates for a newly created
  `"epidataCS"` compared with the same data prepared in earlier
  versions should be very small (e.g., for `data("imdepifit")`
  the mean relative difference of coefficients is 3.7e-08, while the
  `logLik()` is `all.equal()`).
  As an alternative, **rgeos** can still be chosen to do the polygon
  operations.

- The **surveillance**-internal code now depends on
  R >= 2.15.2 (for `nlminb()` `NA` fix of PR#15052,
  consistent `rownames(model.matrix)` of PR#14992,
  `paste0()`, `parallel::mcmapply()`). However, the
  required recent version of **spatstat** (1.34-0, for
  **polyclip**) actually needs R >= 3.0.2, which therefore also
  applies to **surveillance**.

## New Features

- Functions `unionSpatialPolygons()` and
  `intersectPolyCircle()` are now exported. Both are wrappers
  around functionality from different packages supporting polygon
  operations: for determining the union of all subpolygons of a
  `"SpatialPolygons"` object, and the intersection of a polygonal
  and a circular domain, respectively.

- `discpoly()` moved back from [**polyCub**](https://CRAN.R-project.org/package=polyCub)
  to **surveillance**.

## Minor Changes

- **surveillance** now Depends on [**polyCub**](https://CRAN.R-project.org/package=polyCub) (>= 0.4-0)
  and not only Imports it (which avoids `::`-references in
  .GlobalEnv-made functions).

- Nicer default axis labels for `iafplot()`.

- For `twinstim()`, the default is now to `trace`
  every iteration instead of every fifth only.

- Slightly changed default arguments for `plot.epidata()`:
  `lwd` (1->2), `rug.opts` (`col` is set according to
  `which.rug`)

- `twinstim()` saves the vector of `fixed`
  coefficients as part of the returned `optim.args` component,
  such that these will again be held fixed upon `update()`.

- The `plot`-method for `hhh4()`-fits allows for
  region selection by name.


# surveillance 1.6-0 (2013-09-03)

## Synopsis

- The `polyCub`-methods for cubature over polygonal domains
  have been moved to the new dedicated package [**polyCub**](https://CRAN.R-project.org/package=polyCub),
  since they are of a rather general use. The `discpoly()` function
  has also been moved to that package.

- As a replacement for the license-restricted **gpclib** package,
  the **rgeos** package is now used by default
  (`surveillance.options(gpclib=FALSE)`) in generating
  `"epidataCS"` (polygon intersections, slightly slower).
  Therefore, when installing **surveillance** version 1.6-0, the
  system requirements for [**rgeos**](https://CRAN.R-project.org/package=rgeos) have to be met, i.e., GEOS
  must be available on the system. On Linux variants this means
  installing 'libgeos' ('libgeos-dev').

- The improved Farrington method described in Noufaily et al. (2012) is now available as function `farringtonFlexible()`.

- New handling of reference dates in `algo.farrington()` for
  `"sts"` objects with `epochAsDate=TRUE`. Instead of always
  going back in time to the next Date in the `"epoch"` slot, the
  function now determines the *closest* Date. Note that this
  might lead to slightly different results for the upperbound
  compared to previously. Furthermore, the functionality is only
  tested for weekly data (monthly data are experimental). The same
  functionality applies to `farringtonFlexible()`.

- To make the different retrospective modelling frameworks of
  the **surveillance** package jointly applicable, it is now possible
  to convert (aggregate) `"epidataCS"`
  (continuous-time continuous-space data) into an `"sts"` object
  (multivariate time series of counts) by the new function
  `epidataCS2sts`.

- Simulation from `hhh4` models has
  been re-implemented, which fixes a bug and makes it more flexible
  and compatible with a wider class of models.

- The `map`-slot of the `"sts"` class now requires
  `"SpatialPolygons"` (only) instead of
  `"SpatialPolygonsDataFrame"`.

- Re-implementation of `oneStepAhead()` for
  `hhh4`-models with a bug fix, some speed-up and more options.

- Slight speed-up for `hhh4()` fits,
  e.g., by more use of `.rowSums()` and `.colSums()`.

- Crucial speed-up for `twinstim()` fits by more efficient
  code: `mapply`, dropped clumsy `for`-loop in
  `fisherinfo`, new argument `cores` for parallel
  computing via forking (not available on Windows).

## New Features

- Using `tiaf.exponential()` in a `twinstim()` now works
  with `nTypes=1` for multi-type data.

- A legend can be added automatically in `iafplot()`.

- The `untie` methods are now able to produce jittered points
  with a required minimum separation (`minsep`).

- `simulate.ah4` gained a `simplify` argument.

- New `update`-method for fitted `hhh4`-models
  (class `"ah4"`).

- `oneStepAhead()` has more options: specify time range
  (not only start), choose type of start values, `verbose`
  argument.

- `pit()` allows for a list of predictive distributions
  (`pdistr`), one for each observation `x`.

- New spatial auxiliary function `polyAtBorder()`
  indicating polygons at the border (for a `"SpatialPolygons"`
  object).

- `animate.epidataCS()` allows for a `main` title and
  can show a progress bar.

## Minor Changes

- Changed parametrization of `zetaweights()` and completed
  its documentation (now no longer marked as experimental).

- `twinstim(...)$converged` is `TRUE` if
  the optimization routine converged (as before) but contains
  the failure message otherwise.

- Increased default `maxit` for the Nelder-Mead optimizer
  in `hhh4` from 50 to 300, and removed default artificial lower
  bound (-20) on intercepts of epidemic components.

- Renamed returned list from `oneStepAhead`
  (mean->pred, x->observed, params->coefficients,
  variances->Sigma.orig) for consistency, and
  `oneStepAhead()$psi` is only non-`NULL` if we have a
  NegBin model.

- Argument order of `pit()` has changed, which is also
  faster now and got additional arguments `relative` and
  `plot`.

- `twinstim(...)$runtime` now contains the complete
  information from `proc.time()`.

## Bug Fixes

- Fixed a bug in function
  `refvalIdxByDate()` which produced empty reference values
  (i.e. `NA`s) in case the Date entries of `epoch` were not
  mondays. Note: The function works by subtracting `1:b` years from the
  date of the range value and then takes the span `-w:w` around this
  value. For each value in this set it is determined whether the
  closest date in the epoch slot is obtained by going forward or
  backward. Note that this behaviour is now slightly changed compared
  to previously, where we *always* went back in time.

- `algo.farrington()`: Reference values too far back in time
  and hence not being in the `"epoch"` slot of the `"sts"`
  object are now ignored (previously the resulting `NA`s caused the
  function to halt). A warning is displayed in this case.

- `hhh4`: The entry *(5,6)* of the marginal Fisher
  information matrix in models with random intercepts in all three
  components was incorrect.
  If `nlminb` was used as optimizer for the variance parameters
  (using the negative marginal Fisher information as Hessian), this
  could have caused false convergence (with warning) or minimally
  biased convergence (without warning).
  As a consequence, the `"Sigma.cov"` component of the
  `hhh4()` result, which is the inverse of the marginal Fisher
  information matrix at the MLE, was also wrong.

- `untie.matrix()` could have produced jittering greater than
  the specified `amount`.

- `hhh4`: if there are no random intercepts, the
  redundant `updateVariance` steps are no longer evaluated.

- `update.twinstim()` did not work with
  `optim.args=..1` (e.g., if updating a list of models with lapply).
  Furthermore, if adding the `model` component only, the
  `control.siaf` and `optim.args` components were lost.

- `earsC` should now also work with multivariate
  `sts` time-series objects.

- The last week in `data(fluBYBW)` (row index 417) has been
  removed. It corresponded to week 1 in year 2009 and was wrong
  (an artifact, filled with zero counts only).
  Furthermore, the regions in `@map` are now ordered the same as
  in `@observed`.

- Fixed start value of the overdispersion parameter in
  `oneStepAhead` (must be on internal log-scale, not
  reparametrized as returned by `coef()` by default).

- When subsetting `"sts"` objects in time, `@start` was
  updated but not `@epoch`.

- `pit` gave `NA` results if any `x[-1]==0`.

- The returned `optim.args$par` vector in `twinstim()`
  was missing any fixed parameters.

- `hhh4()` did not work with time-varying neighbourhood
  weights due to an error in the internal `checkWeightsArray()`
  function.


# surveillance 1.5-4 (2013-04-21)

- Fixed obsolete `.path.package()` calls.

- Small corrections in the documentation.

- `update.twinstim()` performs better in preserving
  the original initial values of the parameters.

- New pre-defined spatial interaction function
  `siaf.powerlawL()`, which implements a _L_agged power-law
  kernel, i.e. accounts for uniform short-range dispersal.


# surveillance 1.5-2 (2013-03-15)

## New Features

- New method for outbreak detection: `earsC`
  (CUSUM-method described in the CDC Early Aberration Reporting
  System, see Hutwagner et al, 2003).

- Yet another p-value formatting function `formatPval()`
  is now also part of the **surveillance** package.

- `polyCub.SV()` now also accepts objects of classes
  `"Polygon"` and `"Polygons"` for convenience.

## New Features for `twinstim()`

- New spatial interaction function `siaf.powerlaw()`,
  a re-parametrization of the now-deprecated `siaf.lomax()`.

- The temporal `plot`-method for class `"epidataCS"`
  now understands the `add` parameter to add the histogram to an
  existing plot window, and auto-transforms the `t0.Date`
  argument using `as.Date()` if necessary.

- `nobs()` methods for classes `"epidataCS"` and
  `"twinstim"`.

- New argument `verbose` for `twinstim()` which, if
  set to `FALSE`, disables the printing of information messages
  during execution.

- New argument `start` for `twinstim()`, where (some)
  initial parameter values may be provided, which overwrite those in
  `optim.args$par`, which is no longer required (as a naive
  default, a crude estimate for the endemic intercept and zeroes for
  the other parameters are used).

- Implemented a wrapper `stepComponent()` for `step()`
  to perform algorithmic component-specific model selection in
  `"twinstim"` models. This also required the implementation of
  suitable `terms()` and `extractAIC()` methods. The single-step
  methods `add1()` and `drop1()` are also available.

- The `update.twinstim()` method now by default uses the
  parameter estimates from the previous model as initial values for
  the new fit (new argument `use.estimates = TRUE`).

- `as.epidataCS()` checks for consistency of the area of
  `W` and the (now really obligatory) area column in
  `stgrid`.

- `simulate.twinstim()` now by default uses the previous
  `nCircle2Poly` from the `data` argument.

- `direction` argument for `untie.epidataCS()`.

- The `toLatex`-method for `"summary.twinstim"` got
  different defaults and a new argument `eps.Pvalue`.

- New `xtable`-method for `"summary.twinstim"` for
  printing the covariate effects as risk ratios (with CI's and p-values).

## New Features for `hhh4()`

- New argument `hide0s` in the `plot`-method for class
  `"ah4"`.

- New argument `timevar` for `addSeason2formula()`,
  which now also works for long formulae.


# surveillance 1.5-1 (2012-12-14)

- The **surveillance** package is again backward-compatible
  with R version 2.14.0, which is now declared as the minimum
  required version.


# surveillance 1.5-0 (2012-12-12)

## Package Infrastructure

- As requested by the CRAN team, examples now run faster. Some
  are conditioned on the value of the new package option
  `"allExamples"`, which usually defaults to `TRUE` (but is
  set to `FALSE` for CRAN checking, if timings are active).

- Moved some rarely used package dependencies to
  "Suggests:", and also removed some unused packages from there.

- Dropped strict dependence on [**gpclib**](https://CRAN.R-project.org/package=gpclib), which has a
  restricted license, for the **surveillance** package to be clearly
  GPL-2. Generation of `"epidataCS"` objects, which makes use of
  **gpclib**'s polygon intersection capabilities, now requires prior
  explicit acceptance of the **gpclib** license via setting
  `surveillance.options(gpclib = TRUE)`. Otherwise,
  `as.epidataCS()` and `simEpidataCS()` may not be used.

## New Features for `twinstim()`

- Speed-up by memoisation of the `siaf` cubature (using
  the [**memoise**](https://CRAN.R-project.org/package=memoise) package).

- Allow for `nlm`-optimizer (really not recommended).

- Allow for `nlminb`-specific control arguments.

- Use of the expected Fisher information matrix can be disabled
  for `nlminb` optimization.

- Use of the `effRange`-trick can be disabled in
  `siaf.gaussian()` and `siaf.lomax()`. The default
  `effRangeProb` argument for the latter has been changed from
  0.99 to 0.999.

- The `twinstim()` argument `nCub` has been replaced
  by the new `control.siaf` argument list. The old
  `nCub.adaptive` indicator became a feature of the
  `siaf.gaussian()` generator (named `F.adaptive` there) and
  does no longer depend on the `effRange` specification, but uses
  the bandwidth `adapt*sd`, where the `adapt` parameter may be
  specified in the `control.siaf` list in the `twinstim()`
  call. Accordingly, the components `"nCub"` and
  `"nCub.adaptive"` have been removed from the result
  of `twinstim()`, and are replaced by `"control.siaf"`.

- The `"method"` component of the `twinstim()` result
  has been replaced by the whole `"optim.args"`.

- The new `"Deriv"` component of `siaf` specifications
  integrates the "siaf$deriv" function over a polygonal domain.
  `siaf.gaussian()` and `siaf.lomax()` use `polyCub.SV()`
  (with intelligent `alpha` parameters) for this task
  (previously: midpoint-rule with naive bandwidth)

- `scaled` `iafplot()` (default `FALSE`). The
  `ngrid` parameter has been renamed to `xgrid` and is more
  general.

- The `"simulate"` component of `siaf`'s takes an
  argument `ub` (upperbound for distance from the source).

- Numerical integration of spatial interaction functions with an
  `Fcircle` trick is more precise now; this slightly changes
  previous results.

- New S3-generic `untie()` with a method for the
  `"epidataCS"` class (to randomly break tied event times and/or
  locations).

- Renamed `N` argument of `polyCub.SV()` to
  `nGQ`, and `a` to `alpha`, which both have new
  default values.
  The optional polygon rotation proposed by Sommariva &
  Vianello is now also implemented (based on the corresponding MATLAB
  code) and available as the new `rotation` argument.

- The `scale.poly()` method for `"gpc.poly"` is now
  available as `scale.gpc.poly()`. The default return class of
  `discpoly()` was changed from `"gpc.poly"` to
  `"Polygon"`.

- An `intensityplot()`-method is now also implemented for
  `"simEpidataCS"`.

## New Features for `hhh4()`

- Significant speed-up (runs about 6 times faster now, amongst
  others by many code optimizations and by using sparse [**Matrix**](https://CRAN.R-project.org/package=Matrix)
  operations).

- `hhh4()` optimization routines can now be customized for
  the updates of regression and variance parameters separately, which
  for instance enables the use of Nelder-Mead for the variance
  updates, which seems to be more stable/robust as it does
  not depend on the inverse Fisher info and is usually faster.

- The `ranef()` extraction function for `"ah4"` objects
  gained a useful `tomatrix` argument, which re-arranges random
  effects in a unit x effect matrix (also transforming CAR effects
  appropriately).

- Generalized `hhh4()` to also capture parametric
  neighbourhood weights (like a power-law decay). The new function
  `nbOrder()` determines the neighbourhood order matrix
  from a binary adjacency matrix (depends on package [**spdep**](https://CRAN.R-project.org/package=spdep)).

- New argument `check.analyticals` (default `FALSE`)
  mainly for development purposes.

## Bug Fixes

- Fixed sign of observed Fisher information matrix in
  `twinstim`.

- Simulation from the Lomax kernel is now correct (via polar
  coordinates).

- Fixed wrong Fisher information entry for the overdispersion
  parameter in `hhh4`-models.

- Fixed wrong entries in penalized Fisher information wrt the
  combination fixed effects x CAR intercept.

- Fixed indexing bug in penalized Fisher calculation in the case
  of multiple overdispersion parameters and random intercepts.

- Fixed bug in Fisher matrix calculation concerning the relation
  of unit-specific and random effects (did not work previously).

- Improved handling of non-convergent / degenerate solutions during
  `hhh4` optimization. This involves using `ginv()` from
  package [**MASS**](https://CRAN.R-project.org/package=MASS), if the penalized Fisher info is singular.

- Correct labeling of overdispersion parameter in
  `"ah4"`-objects.

- Some control arguments of `hhh4()` have more clear
  defaults.

- The result of `algo.farrington.fitGLM.fast()` now
  additionally inherits from the `"lm"` class to avoid warnings
  from `predict.lm()` about fake object.

- Improved 'NAMESPACE' imports.

- Some additional tiny bug fixes, see the subversion log on
  R-Forge for details.


# surveillance 1.4-2 (2012-08-17)

## Package Infrastructure

- The package is now again compatible with older
  releases of R (< 2.15.0) as intended (by defining `paste0()` in
  the package namespace if it is not found in R **base** at
  installation of the **surveillance** package).

## New Features

- Important new `twinstim()`-feature: fix parameters
  during optimization.

- Useful `update`-method for `"twinstim"`-objects.

- New `[[`- and `plot`-methods for
  `"simEpidataCSlist"`-objects.

- `simEpidataCS()` received tiny bug fixes and is now
  able to simulate from epidemic-only models.

- `R0`-method for `"simEpidataCS"`-objects (actually
  a wrapper for `R0.twinstim()`).

- Removed `dimyx` and `eps` arguments from
  `R0.twinstim()`; now uses `nCub` and
  `nCub.adaptive` from the fitted model and applies the same
  (numerical) integration method.

- `animate.epidata` is now compatible with the
  [**animation**](https://CRAN.R-project.org/package=animation) package.

- More thorough documentation of `"twinstim"`-related
  functions *including many examples*.

## Bug Fixes for `twinstim()`

- `nlminb` (instead of `optim`'s `"BFGS"`) is
  now the default optimizer (as already documented).

- The `twinstim`-argument `nCub` can now be omitted when
  using `siaf.constant()` (as documented) and is internally set to
  `NA_real_` in this case. Furthermore, `nCub` and
  `nCub.adaptive` are set to `NULL` if there is
  no epidemic component in the model.

- `toLatex.summary.twinstim` now again works for
  `summary(*, test.iaf=FALSE)`.

- `print`- and `summary`-methods for
  `"epidataCS"` no longer assume that the `BLOCK` index
  starts at 1, which may not be the case when using a subset in
  `simulate.twinstim()`.

- The `"counter"` step function returned by
  `summary.epidataCS()` does no longer produce false
  numbers of infectives (they were lagged by one timepoint).

- `plot.epidataCS()` now resolves ... correctly and
  the argument `colTypes` takes care of a possible
  `subset`.

- `simEpidataCS()` now also works for endemic-only models
  and is synchronised with `twinstim()` regarding the
  way how `siaf` is numerically integrated (including the
  argument `nCub.adaptive`).

- Fixed problem with `simEpidataCS()` related to missing
  'NAMESPACE' imports (and re-exports) of `marks.ppp` and
  `markformat.default` from [**spatstat**](https://CRAN.R-project.org/package=spatstat), which are required
  for `spatstat::runifpoint()` to work, probably because
  **spatstat** currently does not register its S3-methods.

- Improved error handling in `simEpidataCS()`. Removed a
  `browser()`-call and avoid potentially infinite loop.

## Bug Fixes for `twinSIR()`

- The `.allocate` argument of `simEpidata()` has
  now a fail-save default.

- Simulation without endemic `cox()`-terms now works.

## Minor Changes

- Simplified `imdepi` data to monthly instead of weekly
  intervals in `stgrid` for faster examples and reduced package
  size.

- The environment of all predefined interaction functions for
  `twinstim()` is now set to the `.GlobalEnv`. The previous
  behaviour of defining them in the `parent.frame()` could have
  led to huge `save()`'s of `"twinstim"` objects even with
  `model=FALSE`.

- `simulate.twinSIR` only returns a list of epidemics if
  `nsim > 1`.

- `simulate.twinstim` uses `nCub` and
  `nCub.adaptive` from fitted object as defaults.

- Removed the ...-argument from `simEpidataCS()`.

- The coefficients returned by `simEpidataCS()` are now stored
  in a vector rather than a list for compatibility with
  `"twinstim"`-methods.

- Argument `cex.fun` of `intensityplot.twinstim()` now
  defaults to the `sqrt` function (as in `plot.epidataCS()`.


# surveillance 1.4 (2012-07-26)

## Synopsis

- Besides minor bug fixes, additional functionality has entered the package
  and a new attempt is made to finally release a new version on CRAN
  (version 1.3 has not appeared on CRAN), including a proper 'NAMESPACE'.

## New Features

- Support for non-parametric back-projection using the function
  `backprojNP()` which returns an object of the new
  `"stsBP"` class which inherits from `"sts"`.

- Bayesian nowcasting for discrete time count data is implemented in
  the function `nowcast()`.

- Methods for cubature over polygonal domains can now also visualize what
  they do. There is also a new quasi-exact method for cubature of the
  bivariate normal density over polygonal domains. The
  function `polyCub()` is a wrapper for the different
  methods.

- `residuals.twinstim()` and `residuals.twinSIR()`:
  extract the "residual process", see Ogata
  (1988). The residuals of `"twinSIR"` and
  `"twinstim"` models may be checked graphically by the new
  function `checkResidualProcess()`.

## Significant Changes for `"twinstim"`

- Modified arguments of `twinstim()`: new ordering, new
  argument `nCub.adaptive`, removed argument
  `typeSpecificEndemicIntercept` (which is now specified as part of
  the `endemic` formula as `(1|type)`).

- Completely rewrote the `R0`-method (calculate "trimmed" and
  "untrimmed" *R_0* values)

- The "trimmed" `R0` values are now part of the
  result of the model fit, as well as `bbox(W)`. The
  model evaluation environment is now set as attribute of the
  result if `model=TRUE`.

- New predefined spatial kernel: the Lomax power law kernel
  `siaf.lomax()`

- `plot`-methods for `"twinstim"`
  (`intensityplot()` and `iafplot()`)

- `as.epidataCS()` now auto-generates the stop-column if this is missing

- `print`-method for class `"summary.epidataCS"`

- `[`- and subset-method for `"epidataCS"`
  (subsetting `...$events`)

- `plot`-method for `"epidataCS"`

## Minor Changes

- Improved documentation for the new functionalities.

- Updated references.

- `twinSIR`'s `intensityPlot` is now a method of the
  new S3-generic function `intensityplot`.


# surveillance 1.3 (2011-04-25)

## Synopsis

- This is a major release integrating plenty of new code (unfortunately
  not all documented as good as it could be). This includes code
  for the `"twinstim"` and the `"hhh4"` model.
  The `"twinSIR"` class of models has been
  migrated from package **RLadyBug** to **surveillance**.
  It may take a while before this version will become available from CRAN.

## Significant Changes

- Renamed the `"week"` slot of the `"sts"` S4 class to `"epoch"`.
  All saved data objects have accordingly be renamed, but some hassle
  is to be expected if one you have old `"sts"` objects stored in binary
  form. The function `convertSTS()` can be used to
  convert such "old school" `"sts"` objects.

- Removed the functions `algo.cdc()` and `algo.rki()`.

## New Features

- Support for `"twinSIR"` models (with associated
  `"epidata"` objects) as described
  in H&ouml;hle (2009) has been moved from package
  **RLadyBug** to **surveillance**.
  That means continuous-time discrete-space SIR models.

- Support for `"twinstim"` models as described in
  Meyer et al (2012). That means continuous-time
  continuous-space infectious disease models.

- Added functionality for non-parametric back projection
  (`backprojNP()`) and
  now-casting (`nowcast()`) based on `"sts"` objects.


# surveillance 1.2-2

- Replaced the deprecated `getSpPPolygonsLabptSlots()` calls
  by `sp::coordinates()` when plotting the map slot.

- Minor proof-reading of the documentation.

- Added an argument `"extraMSMargs"` to the algo.hmm function.

- Fixed bug in `outbreakP()` when having observations equal to zero
  in the beginning. Here, $\hat{\mu}^{C1}$ in (5) of Frisen et al (2008)
  is zero and hence the log-based summation in the code failed.
  Changed to product as in the original code, which however might be
  less numerically stable.

- Fixed bug in stcd which added one to the calculated index of idxFA and idxCC.
  Thanks to Thais Rotsen Correa for pointing this out.


# surveillance 1.2-1 (2010-06-10)

- Added `algo.outbreakP()` (Frisen & Andersson, 2009) providing a
  semiparametric approach for outbreak detection for Poisson
  distributed variables.

- Added a pure R function for extracting ISO week and year from Date
  objects. This function (isoWeekYear) is only called if `"%G"` and `"%V"`
  format strings are used on Windows (`sessionInfo()[[1]]$os == "mingw32"`)
  as this is not implemented for `"format.Date"` on Windows.
  Thanks to Ashley Ford, University of Warwick, UK for identifying
  this Windows specific bug.

- For `algo.farrington()` a faster fit routine `"algo.farrington.fitGLM.fast"`
  has been provided by Mikko Virtanen, National Institute for Health
  and Welfare, Finland. The new function calls `glm.fit()`
  directly, which gives a doubling of speed for long series. However, if one
  wants to process the fitted model output some of the GLM routines might
  not work on this output. For backwards compatibility the argument
  `control$fitFun = "algo.farrington.fitGLM"` provides the old (and slow)
  behaviour.


# surveillance 1.1-6 (2010-05-25)

- A few minor bug fixes

- Small improvements in the C-implementation of the `twins()`
  function by Daniel Saban&eacute;s Bov&eacute; fixing the segmentation fault
  issue on 64-bit architectures.


# surveillance 1.1-2 (2009-10-15)

- Added the functions categoricalCUSUM and LRCUSUM.runlength
  for the CUSUM monitoring of general categorical time series
  (binomial, beta-binomial, multinomial, ordered response,
  Bradley-Terry models).

- Added the functions pairedbinCUSUM and pairedbinCUSUM.runlength
  implementing the CUSUM monitoring and run-length computations for
  a paired binary outcome as described in Steiner et al. (1999).

- Experimental implementation of the prospective space-time cluster
  detection described in Assuncao and Correa (2009).

- Added a `demo("biosurvbook")` containing the code of an upcoming
  book chapter on how to use the surveillance package. This
  contains the description of ISO date use, negative binomial CUSUM,
  run-length computation, etc. From an applicational point of view
  the methods are illustrated by Danish mortality monitoring.

- Fixed a small bug in algo.cdc found by Marian Talbert Allen
  which resulted in the control$m argument being ignored.

- The constructor of the sts class now uses the argument
  `"epoch"` instead of weeks to make clearer that also
  daily, monthly or other data can be handled.

- Added additional epochAsDate slot to sts class. Modified
  plot functions so they can handle ISO weeks.

- algo.farrington now also computes quantile and median of
  the predictive distribution. Furthermore has the computation
  of reference values been modified so its a) a little bit faster
  and b) it is also able to handle ISO weeks now. The reference values
  for date t0 are calculated as follows:
  For i, i=1,..., b look at date t0 - i*year. From this date on move
  w months/weeks/days to the left and right. In case of weeks:
  For each of these
  determined time points go back in time to the closest Monday

- Renamed the functions obsinyear to epochInYear, which now also
  handles objects of class Date.


# surveillance 1.0-2 (2009-03-06)

- Negative Binomial CUSUM or the more general NegBin likelihood ratio
  detector is now implemented as part of algo.glrnb.
  This includes the back calculation of the required number of cases
  before an alarm.

- Time varying proportion binomial CUSUM.


# surveillance 0.9-10

- Current status: Development version available from
  <http://surveillance.r-forge.r-project.org/>

- Rewriting of the plot.sts.time.one function to use polygons
  instead of lines for the number of observed cases. Due cause
  a number of problems were fixed in the plotting of the legend.
  Plotting routine now also handles binomial data, where the
  number of observed cases y are stored in `"observed"` and the
  denominator data n are stored in `"populationFrac"`.

- Problems with the aggregate function not operating correctly
  for the populationFrac were fixed.

- The `"rogerson"` wrapper function for algo.rogerson was modified so it
  now works better for distribution `"binomial"`. Thus a time varying
  binomial cusum can be run by calling
  `rogerson( x, control(..., distribution="binomial"))`

- An experimental implementation of the twins model documented in
  Held, L., Hofmann, M., H&ouml;hle, M. and Schmid V. (2006). A two-component
  model for counts of infectious diseases, Biostatistics, 7, pp.
  422--437 is now available as algo.twins.


# surveillance 0.9-9 (2008-01-21)

- Fixed a few small problems
  which gave warnings in the CRAN distribution


# surveillance 0.9-8 (2008-01-19)

- The algo_glrpois function now has an additional `"ret"` arguments,
  where one specifies the return type. The arguments of the underlying
  c functions have been changed to include an additional direction and
  return type value arguments.

- added restart argument to the algo.glrpois control object, which
  allows the user to control what happens after the first alarm has been
  generated

- experimental algo.glrnb function is added to the package. All calls to
  algo.glrpois are now just alpha=0 calls to this function. However,
  the underlying C functions differentiate between poisson and negative case
