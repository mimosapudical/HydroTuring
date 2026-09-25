# Wave-celerity feasibility audit

Base commit: `02da04806195a40e6acbc2750c98748d4aac0b94`

This is an exploratory audit, not a registered HydroTuring probe and not a
proposal.  Its purpose is to answer one question before claiming the ROADMAP
item `momentum/wave-celerity-bounds`:

> Can a transient hydraulic test separate a state-dependent Saint-Venant wave
> from the fixed-celerity routing that already passes the merged
> `momentum/routing-lag-consistency` probe?

## Why this audit exists

The merged routing-lag probe (#82) checks whether travel time has a plausible
magnitude and grows with catchment/channel geometry.  During review, the
`flex_lumped` adapter was deliberately made an independent physical must-pass
model.  It maps channel length to lag using a fixed
`FLOOD_WAVE_CELERITY_M_S = 1.0`, and it passes the routing-lag gate.

That means a new celerity probe is only scientifically distinct if it can catch
a router whose lag scales correctly with length but whose propagation speed
does not respond to hydraulic state.

This audit uses that fixed-celerity behaviour as the negative control.

Relevant merged PR:
https://github.com/Flood-Lab/HydroTuring/pull/82

## Experiment

The audit script is `experiments/wave_celerity_audit.py`.

It imports the existing merged finite-volume model
`models/reference_saint_venant/ht_adapter.py` and calls its real
`solve_steady_reach()` and `_advance()` functions.  No second hydraulic
solver is implemented here.

For each seed:

- draw a wide rectangular reach:
  - width 60–120 m
  - slope 3e-4–1.2e-3
  - Manning n 0.030–0.045
- draw a short reach of 20–30 km and a long reach of 65–85 km;
- use three steady operating depths: 0.5, 1.0 and 2.0 m;
- derive the corresponding steady discharge from the exact rectangular Manning
  rating;
- solve each reach to steady state with the existing Saint-Venant model;
- apply a 5% Gaussian upstream-discharge perturbation;
- integrate the actual Saint-Venant solver forward in time;
- measure the positive-response centroid at the model's centre gauge.

The input pulse is identical in the short and long runs, so source timing
cancels.  The observed propagation speed is

```
c_obs = ((L_long - L_short) / 2) / (t_long - t_short)
```

The independent comparison value is the local kinematic-wave speed

```
c_kin = dQ / dA
```

computed from the exact rectangular Manning rating at the base state.  The
Saint-Venant transient solver does not call or invert this rating to advance
the wave.

The cases remain subcritical: maximum Froude number across the 50-seed main
audit was 0.368.

## Main result: 50 seeds, 64 cells

All **50 / 50** seeds showed

```
c_low < c_medium < c_high
```

for the Saint-Venant reference.

Across all 150 operating points:

- observed/theoretical celerity ratio:
  - min: **0.869**
  - median: **0.887**
  - max: **0.913**
- observed high-flow / low-flow celerity ratio:
  - min: **2.373**
  - median: **2.455**
  - max: **2.491**

Observed celerity ranges across the random geometries:

| state | observed c (m/s) | kinematic theory (m/s) |
|---|---:|---:|
| low | 0.438–0.986 | 0.493–1.097 |
| medium | 0.685–1.559 | 0.775–1.721 |
| high | 1.051–2.438 | 1.209–2.671 |

The systematic low bias is expected from numerical diffusion on the 64-cell
Rusanov grid, so resolution was checked separately rather than hidden inside a
wide tolerance.

## Grid sensitivity

At 128 cells on seeds 0–9:

- 10 / 10 preserve `low < medium < high`;
- observed/theory ratio: **0.923–0.958**, median **0.940**;
- median absolute celerity change relative to 64 cells: about **5.5%**.

At 256 cells on seeds 0–2:

- 3 / 3 preserve `low < medium < high`;
- observed/theory ratio: **0.961–0.981**.

The propagation estimate moves toward the independent kinematic limit as the
grid is refined, which is the desired convergence pattern.  The state
dependence does not disappear.

## Pulse-amplitude sensitivity

For seeds 0–2 at 64 cells, repeating the experiment with 2%, 5% and 10%
upstream pulses produced essentially the same separation:

| pulse | observed/theory range | high/low c ratio |
|---|---:|---:|
| 2% | 0.875–0.911 | 2.408–2.469 |
| 5% | 0.876–0.913 | 2.410–2.472 |
| 10% | 0.876–0.915 | 2.412–2.476 |

So the result is not an artefact of one chosen perturbation amplitude.

## Fixed-celerity negative control

The current `flex_lumped` adapter contains:

```python
FLOOD_WAVE_CELERITY_M_S = 1.0
```

and converts supplied channel geometry to travel time using that fixed speed.
That adapter is a must-pass model for the merged routing-lag probe.

Therefore its state response is, by construction,

```
c_low = c_medium = c_high = 1.0 m/s
```

It can satisfy "longer path -> longer lag" while failing
"hydraulic state -> different wave speed".

This is the exact discrimination a new wave-celerity probe would need in order
not to be a unit-converted duplicate of #82.

## Distinction from neighbouring momentum probes

- `routing-lag-consistency` (#82): checks lag magnitude and geometric scaling.
  Fixed 1 m/s routing passes.
- `routing-causality` (#86 proposal): checks sign / anticipation.  A fixed
  positive celerity is causal and therefore is not the failure targeted here.
- `uniform-flow-friction-consistency` (#114): checks the steady friction
  balance.  A model can use a valid steady Manning rating and still route every
  transient at the same speed.
- `froude-regime` (#103): checks whether reported discharge/depth are in an
  admissible flow regime, not whether the transient characteristic speed
  changes with state.

The proposed scientific gap is therefore specifically **transient
state-dependent propagation**.

## What this audit does *not* establish

This audit does not justify a universal rule `c = 5v/3`.  That approximation
is the wide-channel Manning/kinematic limit, not a law for arbitrary dynamic
waves, backwater, tides or rapidly varied flow.

A production probe would need to construct a regime in which the kinematic
limit is the intended reference condition: mild slope, prismatic rectangular
channel, subcritical flow, no adverse downstream boundary, and a small
perturbation around a steady state.

It also does not freeze a production tolerance.  The 64/128/256-cell results
show why a tolerance must be calibrated against resolution and an independent
physical reference, rather than chosen from the analytic formula alone.

## Decision from this audit

**Feasibility result: positive.**

The candidate is now scientifically stronger than the earlier
"paired length + Manning formula" idea:

1. the positive reference is an independently solved transient Saint-Venant
   system;
2. the negative control is a real repository behaviour already accepted by the
   existing routing-lag probe;
3. the new failure mode is not geometric lag, causality, steady friction or
   Froude regime;
4. separation survives 50 random geometries, grid refinement and pulse-amplitude
   variation.

The next repository-facing step should therefore be a focused proposal issue
for `momentum/wave-celerity-bounds`, **not** implementation yet.  The proposal
should explicitly ask the maintainer to confirm the applicability boundary:
whether fixed-celerity routing is intended to fail this probe whenever a model
claims geometry-aware hydraulic routing, or whether the probe should be limited
to models that declare a Manning/kinematic routing capability.

Until that scope decision is accepted, no production criterion, registry entry,
docs or benchmark archive should be added.
