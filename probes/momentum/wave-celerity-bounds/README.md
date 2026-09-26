# momentum/wave-celerity-bounds

This probe implements accepted proposal #148. It asks a question deliberately
left open by `momentum/routing-lag-consistency`: does transient propagation
speed respond to hydraulic state, rather than merely staying causal and getting
slower in longer reaches?

## Why this is a separate momentum probe

The neighbouring momentum checks constrain different projections of the same
hydrograph, but none identifies the transient channel-wave speed itself.

| Check | What it constrains | What remains invisible |
| --- | --- | --- |
| `momentum/routing-conservation` | routed water is not created or destroyed | timing and propagation speed |
| `momentum/routing-lag-consistency` | gross rainfall-to-runoff lag grows on a geomorphic travel-time scale | runoff generation, storage and channel travel are mixed into one lag; a fixed positive celerity can pass |
| `momentum/stage-discharge-monotonic` | the reported rating is monotone and any hysteresis has the physical sign | how fast a transient moves down the reach |
| `momentum/uniform-flow-friction-consistency` | steady discharge and stage satisfy the declared Manning balance | the characteristic speed of a transient about that steady state |
| `momentum/routing-causality` (accepted proposal) | routed flow does not anticipate its driver | any positive but physically wrong propagation speed |
| `momentum/froude-regime` (accepted proposal) | the reported hydraulic state stays in an admissible regime | state-dependent wave speed inside that regime |

The physical gap is also one order deeper than a new state constraint. The
neighbouring steady probes test values of the hydraulic state or an equilibrium
constitutive relation. #148 tests the **local dynamic derivative** of that
relation: the characteristic speed seen by a small perturbation. A model can
sit on a Manning-consistent steady state and still use an unrelated transient
router. In that sense the probe asks whether the model's transient Jacobian is
consistent with the physics that defines its equilibrium state, not merely
whether the equilibrium itself looks plausible.

The identification step is the contribution as much as the final inequality.
A single rainfall-to-runoff lag cannot separate channel travel from the common
runoff-generation and storage clock. In the controlled paired experiment, write

```
t_response(L, Q) = t_common(Q) + t_route(L, Q).
```

The short and long variants have byte-identical forcing and identical static
attributes except `reach_length_m`. Their difference therefore removes the
common term,

```
Delta t(Q) = t_response(L_long, Q) - t_response(L_short, Q),
```

so `Delta x / Delta t(Q)` identifies the propagation speed associated with
the additional reach length using only the externally reported discharge.
Repeating the same identification at low, medium and high base flow then asks
a question no existing probe asks: whether that isolated propagation speed
changes with hydraulic state as the declared reach physics requires.

A fixed-celerity router is the concrete blind spot. It can conserve mass, be
strictly causal, delay a longer reach more than a shorter one, remain
subcritical, and even have a Manning-consistent steady rating while still using
the same transient speed at every discharge. #148 is designed to reject that
specific construction.

## Hydraulic theory behind the estimator

This probe is a black-box system-identification experiment for the
one-dimensional flood-wave problem, not another empirical lag rule.

With no lateral inflow, continuity is

```
dA/dt + dQ/dx = 0.
```

In the kinematic limit the reach has an equilibrium rating `Q = Q(A)`. For a
small perturbation `a` about a steady state `A0`,

```
A = A0 + a
Q(A) = Q(A0) + (dQ/dA)|A0 a + O(a^2),
```

so the linearized continuity equation is

```
da/dt + c(A0) da/dx = 0,
c(A0) = (dQ/dA)|A0.
```

Thus the Kleitz--Seddon speed `dQ/dA` is the characteristic speed of the
small flood-wave perturbation in this controlled regime; it is not the mean
water velocity and it is not the gross rainfall-to-runoff lag. Mishra & Singh
(2001, https://doi.org/10.1080/02626660109492830) connect the Seddon formula to
the linearized Saint-Venant solution. Ponce & Simons
(1977, https://doi.org/10.1061/JYCEAJ.0004892) show that the full shallow-water
problem contains distinct gravity, dynamic and kinematic wave bands, which is
why this probe deliberately stays in the small, friction-dominated kinematic
regime.

For a wide rectangular Manning reach, `Q proportional to A^(5/3)`, hence

```
c = dQ/dA = (5/3) Q/A = (5/3) u.
```

The production criterion uses the exact rectangular hydraulic radius rather
than this wide-channel approximation.

### Why the response centroid is the measured clock

Keeping the next-order pressure effect gives the linear convection--diffusion
(or diffusion-wave) form

```
dq'/dt + c dq'/dx = D d2q'/dx2.
```

River-routing literature treats `c` and hydraulic diffusivity `D` as the two
physical propagation parameters; the Hayami solution is the classical
constant-parameter analytical case (Moussa 1996,
https://doi.org/10.1002/(SICI)1099-1085(199609)10:9%3C1209::AID-HYP380%3E3.0.CO;2-2).
A modern statement of the same diffusion-wave model gives
`c = dQ/dA` and, in the usual low-inertia approximation,
`D approximately Q/(2 B S)`
(https://doi.org/10.1029/2023WR034692).

For a reach of length `L`, the Hayami impulse kernel is

```
g_L(t) =
  L / (2 sqrt(pi D t^3))
  exp(-(L - c t)^2 / (4 D t)),    t > 0.
```

Its first two temporal moments are

```
E[T_L]   = L / c,
Var[T_L] = 2 D L / c^3.
```

This gives the paired experiment a stronger interpretation than "subtract two
lags". If the common runoff-generation/storage response is `f(t)`, the outlet
response is the convolution `f * g_L`. Temporal cumulants add under
convolution, so

```
centroid(out_L) = centroid(f) + L/c.
```

Therefore

```
centroid(out_long) - centroid(out_short)
    = (L_long - L_short) / c,
```

and the unknown common upstream timing cancels. Diffusion can broaden and
attenuate the hydrograph without moving this first-moment identity. This is why
the probe measures a response centroid rather than a peak index. The second
cumulant also gives a natural future diagnostic,

```
Var(out_long) - Var(out_short)
    = 2 D (L_long - L_short) / c^3,
```

so the same paired design can in principle separate translation (`c`) from
dispersion (`D`) without reading model internals. #148 keeps the verdict on
celerity only.

Temporal moments and cumulants are established tools in flood-routing theory,
not a new law introduced by this probe. Åkesson et al. (2015,
https://doi.org/10.1002/2014WR016279) derived central temporal moments from a
kinematic-diffusive wave description to study stage-dependent hydraulic
response in stream networks, and Romanowicz & Doroszkiewicz (2019,
https://doi.org/10.26491/MHWM/95023) review the use of impulse-response
cumulants for linearized Saint-Venant routing. The contribution here is the
controlled **counterfactual pairing**: evaluate the same black-box model under
the same forcing and hydraulic state at two reach lengths, then difference
those moments so the unknown common land/storage response cancels before the
hydraulic derivative is tested across states.

This interpretation also has an empirical analogue. Allen et al. (2018,
https://doi.org/10.1029/2018GL077914) estimated river-wave celerity from paired
upstream/downstream gauges by dividing reach distance by an observed hydrograph
lag. Meyer et al. (2019, https://doi.org/10.1080/02626667.2018.1557336) found
that celerity--discharge relations vary strongly across real rivers and can
reverse after floodplain activation. That is exactly why #148 restricts itself
to an in-bank prismatic reach instead of asserting monotonic celerity for
arbitrary natural-river states.

## Experiment

For each seed, the generator draws one mild prismatic rectangular reach
geometry and creates three positive base discharges: 8, 16 and 32 m3/s. Only two model runs are needed per seed: a 4 km `short` reach and a 20 km
`long` reach with byte-identical forcing. Inside each run, low, medium and high
hydraulic states occupy successive eight-day blocks. Each block gets four days
to settle, then a six-hour +5% pulse and more than three days of response tail.

Both runs are read at the reach outlet, so the paired propagation distance is

```
Delta x = L_long - L_short.
```

The criterion removes the pre-event discharge, takes the centroid of the
positive transient response, and estimates

```
c_obs = Delta x / (t_long - t_short).
```

For each hydraulic state, differencing the paired reaches cancels the common
forcing and runoff-generation clock. A 72-hour local response window ends
before the next state transition, so the following plateau cannot pull the
previous pulse centroid downstream. The generated record is hourly so timing
quantisation is smaller than the short/long travel-time difference.

## Scope and identification assumptions

The paired estimator is intentionally narrower than a generic flood-wave test.
Its interpretation is valid only when changing `reach_length_m` changes the
propagation distance without changing the hydraulic experiment itself. The
criterion therefore enforces, rather than merely documents, the two conditions
on which the cancellation rests:

- short and long runs carry byte-identical forcing and identical static
  attributes except `reach_length_m`;
- their pre-pulse base discharges agree within one percent.

A model/domain in which changing the declared reach length also changes the
settled base flow has not answered the paired experiment. That is why such a
configuration is not rescued by widening the celerity tolerance. The contract
can verify that a model consumes `reach_length_m`, but it cannot prove where
inside that model the value is used. If reach length is allowed to alter runoff
generation itself, the paired estimate is a net length-dependent response and
must not be interpreted as a pure channel celerity; the intended applicability
is to models that use this field as reach-routing geometry.

The physical comparison is also deliberately a kinematic-limit statement: a
mild prismatic rectangular reach, positive slope and discharge, no adverse
downstream boundary or backwater, and a small perturbation about a settled
state. Passing does not claim that `c = dQ/dA` describes tidal, strongly
backwatered, rapidly varied or supercritical flow. Those regimes require
different probes and boundary information.

## Independent expectation

For a rectangular section,

```
Q = (1/n) A R_h^(2/3) sqrt(S)
R_h = A / P
P = b + 2h
A = b h
```

and differentiating this declared rating gives

```
c_kin = dQ/dA
      = Q [5/(3A) - 4/(3 b P)].
```

The criterion obtains the base discharge from the model output, inverts the
monotone Manning relation for depth, and evaluates this derivative. It does not
reuse a propagation time from the reference model.

The three requirements are:

1. `c_obs > 0`;
2. `abs(c_obs - c_kin) / c_kin <= 0.05`;
3. `c_low < c_medium < c_high`.

The five-percent physical allowance is fixed before gate evaluation. Proposal
#148's feasibility audit found the existing 64-cell Rusanov reference biased
low by about 9--13%, 128 cells by about 4--8%, and tested 256-cell cases by
about 2--4%. The production transient path therefore uses 256 cells rather than
widening the physical tolerance to cover coarse-grid diffusion.

The production generator keeps the synthetic section deliberately wide
(60--90 m for 8--32 m3/s) and uses a 4 km / 20 km pair. Because the reported
discharge is read at the reach outlet, the paired propagation distance is the
full 16 km length difference. The resulting short/long centroid delay remains
resolved at the one-hour output step even at high flow.

## Baselines and applicability

| Model | #148 status | Reason |
| --- | --- | --- |
| `reference_bucket` | N/A | no geometry-aware reach-wave process |
| `flex_lumped` | N/A | existing lag uses catchment channel-length fields and fixed 1 m/s celerity, not this reach contract |
| `flex_topo` | N/A | does not consume this reach geometry |
| `sacsma_snow17` | N/A | does not consume this reach geometry |
| `reference_saint_venant` | gate must-pass | independent finite-volume dynamic-wave solve |
| `reference_fixed_celerity` | gate must-fail | deliberately fixed 1 m/s propagation |


`reference_saint_venant` is the CI must-pass model. Its daily steady path is
unchanged; for sub-daily forcing it advances the same continuity and momentum
equations continuously through each output interval using a 256-cell grid.
The criterion does not call the solver's fluxes or internal wave speeds.

`reference_fixed_celerity` is the deliberately broken control. It transports
the prescribed discharge downstream at exactly 1 m/s for every hydraulic
state. It remains positive and length-dependent, but cannot satisfy the
low/medium/high state response.

The four general physical baselines are intentionally not forced through this
probe. `reference_bucket`, `flex_topo` and `sacsma_snow17` do not declare
consumption of reach length; `flex_lumped` consumes catchment channel-length
fields used by the separate routing-lag probe, not this probe's
`reach_length_m`. Under the harness contract they are therefore
`N/A (INCOMPATIBLE)`, not FAIL.

`reference_saint_venant` is the trusted in-repository numerical reference,
not the independent submitted-model evidence required by
`docs/writing-a-probe.md`. Wflow was audited first but rejected for this role.
Its native river `q_av` was tested separately from total catchment runoff, with
the declared length, width, slope and Manning roughness wired into Wflow's
existing river kinematic wave. The short and long cases nevertheless settled
to different base river discharges before the celerity criterion could be
evaluated. In this one-cell land/river/outlet schematisation, changing reach
geometry changes the hydraulic base state as well as propagation distance, so
it is not the otherwise-identical paired experiment this probe requires. No
Wflow PASS row is claimed.

LISFLOOD was also audited as an independent candidate and rejected
rather than tuned to pass. Wiring the declared geometry into its native
kinematic-wave channel produced, on gate seed 397273707, observed paired
celerities of about 12.9, 43.7 and 63.1 m/s versus Manning expectations of
about 0.66, 0.78 and 1.00 m/s. Refining LISFLOOD's channel routing sub-step
from 3600 s to 300 s changed those values only to about 12.8, 43.4 and
62.7 m/s. The mismatch is therefore not the probe's five-percent allowance
or hourly timing quantisation. No LISFLOOD PASS row is claimed or archived.

The independent submitted-physical-model evidence required by
`docs/writing-a-probe.md` remains the outstanding merge prerequisite. No
unimplemented proposal is counted as that evidence. Accepted model proposal
#132, mizuRoute, is scientifically relevant because it is a routing-only
physical model with native hourly reach discharge and explicit kinematic-wave
schemes, but it consumes lateral runoff rather than the rainfall forcing used
by this paired experiment and its separate reach-routing infrastructure is
still under review. #148 therefore does not depend on, claim a result from, or
take ownership of that model contribution. A future submitted model counts
only after it can run this experiment honestly and archive a passing row.

## Reproduction

```sh
ht validate
ht gate --probe momentum/wave-celerity-bounds
pytest -q
```

No generated forcing or model output is committed.
