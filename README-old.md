# Smoothing priors

This package implements smoothing priors and related utilities for dealing with skyline objects in BEAST2.

## What is a skyline?

A skyline is a piecewise constant function on $n$ intervals, $\{f_0, f_1, ..., f_n \}$, defined on $n+1$ time intervals, $\{t_0,t_1, ... , t_{n+1} \}$, such that:

\[
  f(t) = f_i \iff t \in [t_i,t_{i+1})
\]

Thus, a skyline is a type of time-series. There is no requirement for intervals to be the same length. Skylines can run forward or backward in time (i.e. $t_0$ or $t_{n+1}$ may be defined as the present).

Usually a skyline is defined between the present and the tMRCA of a tree or the origin of an epidemic. In birth-death models, the value of the skyline function represents rates of the birth-death model through time or probabilities. In coalescent models the value of the skyline function represents the effective population size ($N_e$).

Within BEAST2 a skyline is defined by two $n$-dimensional parameters, representing the skyline-values and the shift times, e.g.:

```xml
  <parameter id="reproductiveNumber" dimension="5">2.0</parameter>
  <parameter name="intervalTimes" id="intervalTimes" value="0 0.2 0.4 0.6 0.8"/>  
```

The $(n+1)^{th}$ time is always implicit, depending on the model specification used. If relative times are used it would be 1.0 in the example above. If times are absolute it would probably be the origin or tMRCA in the example above (but this depends on the exact model specifition). Note that if the model has a `reverseTimes` parameter then the implicit time may be 0, depending on the value of `reverseTimes` for the parameter represented by the skyline function.

## What is a smoothing prior?

The specification of a skyline does not by itself impose any dependency between the value of the skyline in different time intervals. Thus, each $f_i$ is assumed to be i.i.d. (This is the default behaviour in BEAST2, when simply applying a single prior function on the whole parameter vector. However, other model components could be adding an implicit dependence between successive entries in the vector).

A smoothing prior makes explicit a dependency between successive entries in the skyline function. Thus, we can have a prior on the first (or last) entry of the vector and have some relationship between each successive entry. This makes sense for parameters where we expect it to change smoothly through time, such as the effective reproductive number ($R_e$) or effective population size ($N_e$).

BEAST1.8 and BEAST2 implement the Markov chain distribution smoothing prior in `beast.math.distributions.MarkovChainDistribution`, where a Gamma/LogNormal distribution is used to model the dependency between successive entries. BEAST1.8 also implements a time-aware GMRF smoothing-prior, however it is not implemented in a general way, but only as part of skyride/skygrid models.

**This package aims to implement a variety smoothing priors that can be used in a general way for any skyline/time-series!**

## Why slice trees?

TreeSlicer can be used to create a time vector for a skyline based on events in a tree, or calendar dates, without having to know _a priori_ what the date of most recent/oldest sample or tMRCA of the tree is. This can be very useful for specifying skylines without having to laboriously write down the times of events.

TreeSlicer always returns times from the most recent sample in a tree. Thus, time in the present is always 0.

More detailed help will be added later.
