#!/usr/bin/env python3
"""Generate deterministic Aróstica pressure and active-tension tables."""

from __future__ import annotations

import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def activation(time, systole, diastole, gamma=0.005, minimum=-30.0, maximum=5.0):
    rise = 0.5 * (1.0 + math.tanh((time - systole) / gamma))
    fall = 0.5 * (1.0 - math.tanh((time - diastole) / gamma))
    value = rise * fall
    return maximum * value + minimum * (1.0 - value)


def rk4(value, time, step, derivative):
    k1 = derivative(time, value)
    k2 = derivative(time + step / 2.0, value + step * k1 / 2.0)
    k3 = derivative(time + step / 2.0, value + step * k2 / 2.0)
    k4 = derivative(time + step, value + step * k3)
    return value + step * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0


def integrate(derivative, step=1e-5, sample=1e-3):
    result = []
    value = 0.0
    time = 0.0
    next_sample = 0.0
    while time < 1.0 - 1e-14:
        while next_sample <= time + 1e-14:
            result.append((next_sample, value))
            next_sample += sample
        local_step = min(step, 1.0 - time)
        value = rk4(value, time, local_step, derivative)
        time += local_step
    result.append((1.0, value))
    return result


def active_tension(systole, diastole, sigma0=1.5e5):
    def derivative(time, value):
        signal = activation(time, systole, diastole)
        return -abs(signal) * value + sigma0 * max(signal, 0.0)
    return integrate(derivative)


def chamber_pressure(apre, amid, sigpre, sigmid, systole=0.17, diastole=0.484):
    def derivative(time, value):
        signal = activation(time, systole, diastole)
        preload = 0.5 * (1.0 - math.tanh((time - diastole) / 0.005))
        decay = signal + apre * preload + amid
        return -abs(decay) * value + sigmid * max(decay, 0.0) + sigpre * abs(preload)
    return integrate(derivative)


def write_table(path, values):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as stream:
        stream.write("(\n")
        for time, value in values:
            stream.write(f"({time:.6f} {value:.12g})\n")
        stream.write(")\n")


def main():
    tables = {
        "tau_monoventricle.dat": active_tension(0.16, 0.484),
        "tau_biventricle.dat": active_tension(0.163, 0.5),
        "pressure_caseB.dat": chamber_pressure(5.0, 1.0, 7000.0, 16000.0),
        "pressure_LV.dat": chamber_pressure(5.0, 15.0, 12000.0, 16000.0),
        "pressure_RV.dat": chamber_pressure(1.0, 10.0, 3000.0, 4000.0),
    }
    maxima = {name: max(value for _, value in values) for name, values in tables.items()}
    for case in ("monoventricle", "biventricle"):
        output = ROOT / "cases" / case / "constant" / "loadCurves"
        for name, values in tables.items():
            write_table(output / name, values)
        (output / "maxima.json").write_text(json.dumps(maxima, indent=2) + "\n")
    print(json.dumps(maxima, indent=2))


if __name__ == "__main__":
    main()

