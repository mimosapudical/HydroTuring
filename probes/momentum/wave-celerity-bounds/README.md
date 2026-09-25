# momentum/wave-celerity-bounds

A routing model can be causal, use the supplied reach length, and make a longer
reach respond later while still moving every flood wave at one fixed speed.
This probe asks the missing momentum question: **does transient propagation
speed change with the hydraulic state of the reach?**

## Experiment

For every seed the generator makes one mild, prismatic rectangular reach and
three operating states: low, medium and high flow. Each state is run twice,
with a 10 km and a 50 km reach. Weather, catchment attributes, cross-section,
slope, roughness and the small two-hour perturbation are byte-identical inside
each short/long pair; only reach_length_m changes.

The gauge is at the centre of the reach. If the response centroids in a state
are t_s and t_l, the common runoff-generation delay cancels and the additional
travel distance is

~~~text
Δx = (L_long - L_short) / 2
c_obs = Δx / (t_long - t_short).
~~~

The response centroid is formed from the positive discharge anomaly above the
six-hour pre-pulse median. Values below 1e-4 of the response peak are discarded
before the centroid is formed so that a floating-point tail cannot move a
timing statistic arbitrarily far.

## Independent hydraulic expectation

In the controlled friction-dominated limit, Manning's rectangular rating is

~~~text
Q = (1/n) A R^(2/3) S^(1/2)
R = A / P
P = w + 2 A/w.
~~~

Differentiating that rating, rather than using the wide-channel shortcut, gives

~~~text
c_kin = dQ/dA
      = Q [5/(3A) - 4/(3 w P)].
~~~

The criterion requires, for every state:

1. t_long > t_short, hence downstream c_obs > 0;
2. abs(c_obs - c_kin) / c_kin <= 0.05;
3. the observed celerities increase materially from low to medium to high,
   with each adjacent gap exceeding 2% of the larger speed.

The third condition is the discrimination that the existing
momentum/routing-lag-consistency probe does not make. A constant positive
celerity can produce a perfectly causal length-dependent lag and still fail
this state-response check.

## Numerical resolution and tolerance

The production reference uses 256 finite-volume cells and hourly output on
10 km / 50 km reaches. The five-percent physical allowance was fixed before
the gate-seed run; it was not widened to admit the reference.

A pre-PR numerical audit using the same Rusanov Saint-Venant update as the
repository reference produced the following ratios on the three deterministic
gate seeds:

| gate seed | low c_obs/c_kin | medium | high |
| ---: | ---: | ---: | ---: |
| 397273707 | 0.9727 | 0.9768 | 0.9857 |
| 904225821 | 0.9746 | 0.9780 | 0.9919 |
| 1411177935 | 0.9750 | 0.9761 | 0.9933 |

The largest discrepancy is about 2.73%, leaving more than two percentage
points between the numerical reference error and the 5% gate. All nine state
comparisons satisfy low < medium < high. The generated states remain in the
mild subcritical regime; the probe deliberately does not claim that the local
kinematic limit applies to backwater, tidal, rapidly varied or supercritical
flow.

reference_saint_venant is numerically independent of the criterion. Its
transient mode starts from an equilibrium obtained by integrating the full
one-dimensional continuity and momentum equations, then advances depth and
unit discharge with the Rusanov finite-volume flux and Manning source term.
It never evaluates dQ/dA or a normal-depth formula to decide the transient
arrival time.

## Controls and applicability

The deliberately broken reference_fixed_celerity converts the same effective
forcing to discharge and transports it to the centre gauge at exactly 1 m/s at
all three operating states. Its steady stage is Manning-consistent, so it
isolates the intended defect: the propagation speed does not respond to state.
It is listed under must_fail_only, so the gate requires exactly the new
criterion to catch it.

The repository's four ordinary physical baselines do not currently implement
this explicit reach-transient contract:

| model | current result on this probe | reason |
| --- | --- | --- |
| reference_bucket | N/A (INCOMPATIBLE) | does not consume reach_length_m |
| flex_lumped | N/A (INCOMPATIBLE) | routing uses main_channel_length_km / centroid_channel_length_km, not reach_length_m |
| flex_topo | N/A (INCOMPATIBLE) | does not consume reach_length_m |
| sacsma_snow17 | N/A (INCOMPATIBLE) | does not consume reach_length_m |

This is intentional rather than relabelling an input just to make a gate run.
In particular, flex_lumped contains a fixed 1 m/s routing assumption and is
useful evidence for the failure mode, but changing its declared contract would
change the meaning of an existing physical baseline.

The acceptance guide has a special rule for a process none of those four
models implements: an exact reference may gate it only when an independent
submitted physical model with that process has a passing archived row. That
independent submitted-model row is therefore a **merge prerequisite**, not
something this README treats as already satisfied. It must be produced by a
real compatible router; no result row is fabricated here.

## Scope

Passing establishes one necessary transient momentum property in a controlled
kinematic limit. It does not establish whole-reach mass closure, arbitrary
Saint-Venant characteristic accuracy, backwater behaviour or flood-wave
causality by itself. Those remain separate questions.

## Sources

- MacDonald et al. (1997), Analytic Benchmark Solutions for Open-Channel
  Flows, Journal of Hydraulic Engineering,
  https://doi.org/10.1061/(ASCE)0733-9429(1997)123:11(1041).
- Delestre et al. (2013), SWASHES: a compilation of shallow water analytic
  solutions for hydraulic and environmental studies,
  https://doi.org/10.1002/fld.1839.
- Kurganov (2018), Finite-volume schemes for shallow-water equations,
  Acta Numerica, https://doi.org/10.1017/S0962492918000028.
- USGS Water-Resources Investigations Report 83-4251, rectangular-section
  Manning definitions, https://doi.org/10.3133/wri834251.
