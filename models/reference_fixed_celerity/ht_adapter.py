#!/usr/bin/env python3
"""Deliberately broken fixed-celerity control for wave-celerity bounds."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


MODEL = {"name": "reference_fixed_celerity", "version": "1.0.0"}
CELERITY_M_S = 1.0
SECONDS_PER_DAY = 86400.0
TIMESTEP_SECONDS = {"PT1H": 3600.0}


def _capacity(depth: float, width: float, slope: float, roughness: float) -> float:
    area = width * depth
    radius = area / (width + 2.0 * depth)
    return area * radius ** (2.0 / 3.0) * math.sqrt(slope) / roughness


def _normal_depth(discharge: float, width: float, slope: float, roughness: float) -> float:
    if discharge <= 0.0:
        raise ValueError("fixed-celerity control requires positive discharge")
    lo, hi = 0.0, 0.25
    while _capacity(hi, width, slope, roughness) < discharge:
        hi *= 2.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if _capacity(mid, width, slope, roughness) < discharge:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _delayed(values: list[float], delay_steps: float) -> list[float]:
    """Fractionally delay a regular series, holding its initial base flow."""
    if delay_steps < 0.0 or not math.isfinite(delay_steps):
        raise ValueError("delay must be finite and non-negative")
    out: list[float] = []
    for i in range(len(values)):
        source = i - delay_steps
        if source <= 0.0:
            out.append(values[0])
            continue
        left = int(math.floor(source))
        if left >= len(values) - 1:
            out.append(values[-1])
            continue
        fraction = source - left
        out.append((1.0 - fraction) * values[left] + fraction * values[left + 1])
    return out


def simulate(forcing: list[dict], static: dict, step_s: float) -> tuple[list[dict], dict]:
    area_km2 = float(static["area_km2"])
    width = float(static["width_m"])
    bed = float(static["bed_elevation_m"])
    slope = float(static["slope"])
    roughness = float(static["manning_n"])
    length = float(static["reach_length_m"])
    shape = str(static["cross_section_shape"]).strip().lower()
    if shape != "rectangular":
        raise ValueError("reference_fixed_celerity requires a rectangular section")
    if min(area_km2, width, slope, roughness, length) <= 0.0:
        raise ValueError("area and reach geometry must be positive")

    upstream: list[float] = []
    for item in forcing:
        effective = max(float(item["pr"]) - float(item["pet"]), 0.0)
        upstream.append(effective * 1.0e-3 * area_km2 * 1.0e6 / SECONDS_PER_DAY)

    # Gauge is at the centre of the reach, exactly as in the criterion.
    delay_s = 0.5 * length / CELERITY_M_S
    routed = _delayed(upstream, delay_s / step_s)
    rows = []
    for item, discharge in zip(forcing, routed):
        depth = _normal_depth(discharge, width, slope, roughness)
        rows.append({
            "time": item["time"],
            "dis": discharge,
            "stage": bed + depth,
        })
    return rows, {
        "method": "deliberately fixed celerity",
        "celerity_m_s": CELERITY_M_S,
        "gauge_distance_m": 0.5 * length,
        "delay_s": delay_s,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--request", required=True)
    args = parser.parse_args()
    request_path = Path(args.request).resolve()
    request = json.loads(request_path.read_text())
    io_dir = request_path.parent
    timestep = str(request.get("timestep", "PT1H"))
    if timestep not in TIMESTEP_SECONDS:
        raise SystemExit(f"unsupported timestep {timestep!r}")

    with open(io_dir / request["input"]["forcing"], newline="") as fh:
        forcing = list(csv.DictReader(fh))
    static = json.loads((io_dir / request["input"]["static"]).read_text())
    rows, routing = simulate(forcing, static, TIMESTEP_SECONDS[timestep])

    output = io_dir / request["output"]["table"]
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["time", "dis", "stage"])
        writer.writeheader()
        writer.writerows(rows)
    (io_dir / request["output"]["run"]).write_text(json.dumps({
        "status": "ok",
        "model": MODEL,
        "n_steps": len(rows),
        "routing": routing,
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
