"""Fission-product decay heat — full ANS/ANSI-5.1 23-group model.

Each group ``i`` obeys ``dh_i/dt = alpha_i * F(t) - lambda_i * h_i``, where
``F(t)`` is the fission rate normalized to rated power (the kinetics variable
``n``). The decay-power fraction is ``sum(h_i) / Q`` with ``Q = 200 MeV/fission``
(ANS normalization convention).

Groups are advanced with the exact exponential solution for a constant ``F``
over the substep, so the integration is unconditionally stable at any ``dt``
(including 600x time compression)::

    h_i(t+dt) = h_i * e^(-lambda_i*dt) + (alpha_i/lambda_i) * F * (1 - e^(-lambda_i*dt))

Decay heat sits at ~6.5% of rated power during sustained operation and, after a
trip (``n -> 0``), decays along the ANS curve — which is the entire point of the
model: it sustains the residual power that the plant must still remove once the
chain reaction has stopped.
"""

from __future__ import annotations

import math

# ANS-5.1 U-235 thermal-fission 23-group fit.
# alpha_i [MeV/(fission*s)], lambda_i [1/s].
_ALPHA = (
    6.5057e-01, 5.1264e-01, 2.4384e-01, 1.3850e-01, 5.5440e-02,
    2.2225e-02, 3.3088e-03, 9.3015e-04, 8.0943e-04, 1.9567e-04,
    3.2535e-05, 7.5595e-06, 2.5232e-06, 4.9948e-07, 1.8531e-07,
    2.6608e-08, 2.2398e-09, 8.1641e-12, 8.7797e-11, 2.5131e-14,
    3.2176e-16, 4.5038e-17, 7.4791e-17,
)
_LAMBDA = (
    2.2138e+01, 5.1587e-01, 1.9594e-01, 1.0314e-01, 3.3656e-02,
    1.1681e-02, 3.5870e-03, 1.3930e-03, 6.2630e-04, 1.8906e-04,
    5.4988e-05, 2.0958e-05, 1.0010e-05, 2.5438e-06, 6.6361e-07,
    1.2290e-07, 2.7213e-08, 4.3714e-09, 7.5780e-10, 2.4786e-10,
    2.2384e-13, 2.4600e-14, 1.5699e-14,
)

# Recoverable energy per fission [MeV], used for normalization.
_Q_FISSION = 200.0

# Equilibrium (infinite-operation) decay-heat fraction at rated power:
#   sum(alpha_i / lambda_i) / Q  ~= 0.0654.
EQUILIBRIUM_FRACTION = sum(a / l for a, l in zip(_ALPHA, _LAMBDA)) / _Q_FISSION


class DecayHeat:
    """23-group ANS-5.1 fission-product decay-heat model.

    Parameters
    ----------
    initial_power:
        Normalized fission power (``1.0`` == rated) the core is assumed to have
        been operating at long enough to reach decay-heat equilibrium. The group
        inventories are seeded to that equilibrium, so a plant that starts at
        steady power already carries its full ~6.5% decay-heat load and shows
        correct residual heat immediately after a trip — rather than building it
        up from zero over the first hours of simulated time.
    """

    __slots__ = ("_h", "fraction")

    def __init__(self, initial_power: float = 1.0) -> None:
        f0 = max(initial_power, 0.0)
        # Equilibrium inventory for a constant fission rate f0: h_i = alpha_i/lambda_i * f0.
        self._h = [a / l * f0 for a, l in zip(_ALPHA, _LAMBDA)]
        self.fraction = sum(self._h) / _Q_FISSION

    def step(self, dt: float, n: float) -> float:
        """Advance by ``dt`` seconds at normalized fission power ``n``.

        Returns the decay-heat fraction of rated power after the step.
        """
        f = n if n > 0.0 else 0.0
        h = self._h
        total = 0.0
        for i in range(23):
            lam = _LAMBDA[i]
            e = math.exp(-lam * dt)
            h[i] = h[i] * e + _ALPHA[i] / lam * f * (1.0 - e)
            total += h[i]
        self.fraction = total / _Q_FISSION
        return self.fraction
