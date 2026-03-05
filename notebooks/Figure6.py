#!/usr/bin/env python3
"""
Monte‑Carlo sensitivity analysis of flawless‑copy decay
with an exponentially increasing step‑error rate and Ek uncertainty,
plotting only the mean trajectories and starting at a_i = 1.
"""

import numpy as np
import matplotlib.pyplot as plt

# ── Nominal parameters ────────────────────────────────────────────────
epsilons        = [0.01, 0.015, 0.02, 0.03, 0.05,
                   0.06, 0.08, 0.10, 0.20, 0.50]      # baseline ε₀ at a_i = 1
k_exp_growth    = 0.02                                # exponential growth constant k
ai_range        = np.arange(1, 121)                   # a_i = 1 … 120  ← changed
N_total         = 6.022e23                            # Avogadro’s number

# ── Ek uncertainty model ──────────────────────────────────────────────
Ek_mu           = 0.0
Ek_sigma        = 0.005                               # 1 σ of Ek
n_sim           = 5000                                # Monte‑Carlo draws per ε₀
rng             = np.random.default_rng(42)

# ── Storage for summary statistics ────────────────────────────────────
N_mean  = np.zeros((len(epsilons), ai_range.size))
N_low   = np.zeros_like(N_mean)                       # 2.5 percentile
N_high  = np.zeros_like(N_mean)                       # 97.5 percentile

# ── Monte‑Carlo loop ─────────────────────────────────────────────────
for j, eps0 in enumerate(epsilons):
    # Exponentially rising baseline profile, now starting at a_i = 1
    eps_base = eps0 * np.exp(k_exp_growth * ai_range)

    # One Ek sample per trajectory
    Ek_samples = rng.normal(Ek_mu, Ek_sigma, n_sim)

    # Combine and clip to [0, 1]
    eps_profile = np.clip(eps_base + Ek_samples[:, None], 0.0, 1.0)

    # Survival probability per step and cumulative product over steps
    surv_prob = 1.0 - eps_profile
    N_sim     = N_total * np.cumprod(surv_prob, axis=1)

    # Aggregate statistics
    N_mean[j] = N_sim.mean(axis=0)
    N_low[j]  = np.percentile(N_sim, 2.5, axis=0)
    N_high[j] = np.percentile(N_sim, 97.5, axis=0)

# ── Plotting ─────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(8, 5))

for j, eps0 in enumerate(epsilons):
    ax.plot(ai_range, N_mean[j], label=f"$\\varepsilon_0$ = {eps0}")

ax.set_xlabel(r"Assembly index $a_i$")
ax.set_ylabel(r"Copy number $N$")
ax.set_xticks(np.arange(0, 121, 5), minor=True)
ax.set_xlim(1, 120)                                   # ensure axis starts at 1
ax.tick_params(axis='x', which='minor', length=3)
ax.legend(title="", ncol=2, frameon=False)

ax.grid(True, linestyle="--", linewidth=0.5)
ax.grid(True, which="minor", linestyle=":", linewidth=0.4, alpha=0.6)

plt.tight_layout()
plt.show()