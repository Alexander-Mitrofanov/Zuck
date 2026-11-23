"""Estimate shotgun base depth from a hash sample of canonical k-mer types.

The spectrum contains one entry per sampled sequence type, grouped by the
number of times that sequence was observed. The exposure is the sum of the
expected error-free, eligible k-mer windows divided by input bases. This
converts a Poisson rate for usable windows into base sequencing depth.
"""

from bisect import bisect_left
from functools import lru_cache
import math


FIRST = 2
LAST = 100
COUNTS = tuple(range(FIRST, LAST + 1))
BACKGROUND_FRACTION = 0.02
BACKGROUND = tuple((1 / n) / sum(1 / j for j in COUNTS) for n in COUNTS)
MODEL_VERSION = "truncated-poisson-1"


def _poisson(rate):
    """Poisson probabilities conditional on observing between 2 and 100 copies."""
    # The common P(X=2) factor cancels, including when the rate is very small.
    probabilities = [1.0]
    for copies in range(FIRST + 1, LAST + 1):
        probabilities.append(probabilities[-1] * rate / copies)
    total = sum(probabilities)
    return [value / total for value in probabilities]


@lru_cache(maxsize=1)
def _mean_table():
    # E[X | 2 <= X <= 100] is strictly increasing in the Poisson rate.
    # A fine geometric grid keeps inversion stable near both truncation limits.
    rates = [0.025]
    while rates[-1] < 160.0:
        rates.append(rates[-1] * 1.02)
    means = [sum(n * p for n, p in zip(COUNTS, _poisson(rate))) for rate in rates]
    return means, rates


def _rate_for_mean(mean):
    means, rates = _mean_table()
    index = bisect_left(means, mean)
    if index == 0:
        return rates[0]
    if index == len(rates):
        return rates[-1]
    fraction = (mean - means[index - 1]) / (means[index] - means[index - 1])
    return rates[index - 1] + fraction * (rates[index] - rates[index - 1])


def _fit(observed, rates, weights):
    """EM on the conditional spectrum with a fixed nonspecific background."""
    for _ in range(120):
        probabilities = [_poisson(rate) for rate in rates]
        signal = [(1 - BACKGROUND_FRACTION) * weight for weight in weights]
        assigned = [0.0] * len(rates)
        assigned_copies = [0.0] * len(rates)
        for index, count in enumerate(observed):
            if count == 0:
                continue
            terms = [w * probs[index] for w, probs in zip(signal, probabilities)]
            total = sum(terms) + BACKGROUND_FRACTION * BACKGROUND[index]
            for component, term in enumerate(terms):
                mass = count * term / total
                assigned[component] += mass
                assigned_copies[component] += COUNTS[index] * mass
        total_assigned = sum(assigned)
        if any(mass <= 0 for mass in assigned):
            break
        new_rates = [_rate_for_mean(copies / mass) for copies, mass in zip(assigned_copies, assigned)]
        new_weights = [mass / total_assigned for mass in assigned]
        change = max(abs(a - b) / max(1.0, a) for a, b in zip(rates, new_rates))
        weight_change = max(abs(a - b) for a, b in zip(weights, new_weights))
        rates, weights = new_rates, new_weights
        if max(change, weight_change) < 1e-5:
            break
    spectrum = [BACKGROUND_FRACTION * value for value in BACKGROUND]
    for rate, weight in zip(rates, weights):
        for i, probability in enumerate(_poisson(rate)):
            spectrum[i] += (1 - BACKGROUND_FRACTION) * weight * probability
    likelihood = sum(p * math.log(q) for p, q in zip(observed, spectrum) if p)
    distance = sum(abs(p - q) for p, q in zip(observed, spectrum)) / 2
    return likelihood, distance, sorted(zip(rates, weights))


def _fit_mixture(observed, single_rate):
    # Separated starting points reduce the risk of both components converging
    # to the same dominant peak in an unevenly sequenced community.
    candidates = [([max(0.025, single_rate / 2), min(160.0, single_rate * 1.8)], [0.5, 0.5])]
    mass = copies = 0.0
    overall_mean = sum(n * p for n, p in zip(COUNTS, observed))
    targets = iter((0.2, 0.5, 0.8))
    target = next(targets)
    for n, p in zip(COUNTS, observed):
        mass += p
        copies += n * p
        if mass >= target and mass < 1 - 1e-9:
            rates = [_rate_for_mean(copies / mass), _rate_for_mean((overall_mean - copies) / (1 - mass))]
            candidates.append((rates, [mass, 1 - mass]))
            while target <= mass:
                target = next(targets, math.inf)
    return max((_fit(observed, rates, weights) for rates, weights in candidates), key=lambda result: result[0])


def predict(scan):
    """Return an estimate for one uniform sample or a mixed/uncertain decision.

    Component weights refer to repeated sampled k-mer types inside the fitted
    interval. They are neither species abundances nor fractions of input reads.
    """
    result = {"model": MODEL_VERSION, "depth_x": None, "regime": "uncertain"}
    histogram = scan.get("histogram")
    exposure = scan.get("correct_window_exposure")
    if (not isinstance(histogram, (list, tuple)) or len(histogram) != 257
            or any(type(n) is not int or n < 0 for n in histogram) or histogram[0] != 0):
        return {**result, "reason": "invalid_histogram"}
    if not isinstance(exposure, (float, int)) or not 0.15 <= exposure <= 1:
        return {**result, "reason": "inadequate_usable_window_exposure"}
    fitted = sum(histogram[FIRST:LAST + 1])
    overflow = sum(histogram[LAST + 1:])
    result.update(fitted_repeated_types=fitted, repeated_types_above_100=overflow,
                  correct_window_exposure=exposure)
    if fitted < 350 or sum(histogram[3:LAST + 1]) < 30:
        return {**result, "reason": "insufficient_repeated_types"}
    if overflow / (fitted + overflow) > 0.10:
        return {**result, "reason": "excess_multiplicity_above_supported_range"}
    observed = [count / fitted for count in histogram[FIRST:LAST + 1]]
    start = _rate_for_mean(sum(n * p for n, p in zip(COUNTS, observed)))
    single_ll, single_tv, ((rate, _),) = _fit(observed, [start], [1.0])
    mixture_ll, mixture_tv, components = _fit_mixture(observed, rate)
    single_depth = rate / exposure
    low_depth, high_depth = [r / exposure for r, _ in components]
    gain = mixture_ll - single_ll
    result.update(single_depth_x=single_depth, single_tv=single_tv, mixture_tv=mixture_tv,
                  mixture_gain_nats_per_type=gain,
                  mixture_components=[{"depth_x": r / exposure, "repeated_type_fraction": w}
                                      for r, w in components])
    mixed = (min(w for _, w in components) >= 0.125 and high_depth / low_depth >= 2.1
             and high_depth - low_depth >= 1.5 and gain >= 0.035)
    if mixed:
        if mixture_tv > 0.20 or low_depth < 0.4 or high_depth > 80:
            return {**result, "reason": "heterogeneous_fit_outside_supported_range"}
        return {**result, "regime": "mixed", "reason": "separated_depth_components"}
    if single_tv > 0.20 or not 0.4 <= single_depth <= 80:
        return {**result, "reason": "uniform_fit_outside_supported_range"}
    regime = "low" if single_depth < 5 else "intermediate" if single_depth < 12 else "high"
    return {**result, "depth_x": single_depth, "regime": regime, "reason": "uniform_depth_fit"}
