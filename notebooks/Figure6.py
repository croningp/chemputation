import numpy as np
import matplotlib.pyplot as plt

# --- Parameters -------------------------------------------------------
epsilons = [0.01, 0.015, 0.02, 0.03, 0.05,
            0.06, 0.08, 0.10, 0.20, 0.50]      # step-error rates
N_total   = 6.022e23                            # Avogadro’s number
ai_range  = np.arange(0, 121)                   # a_i = 0 … 120
# ----------------------------------------------------------------------

fig, ax = plt.subplots(figsize=(8, 5))

for eps in epsilons:
    N_perfect = N_total * (1 - eps) ** ai_range
    ax.plot(ai_range, N_perfect, label=f'ε = {eps}')

# Axis labels and title
ax.set_xlabel('Assembly index $a_i$')
ax.set_ylabel('Perfect copies $N_{perfect}$')
ax.set_title('Flawless copy number vs. assembly index\n'
             r'$N_{perfect}=N_{total}(1-\varepsilon)^{a_i}$'
             '  (N_total = 6.022×10$^{23}$)')

# Minor x-ticks every 5
ax.set_xticks(np.arange(0, 121, 5), minor=True)
ax.tick_params(axis='x', which='minor', length=3)

# Legend overlay, transparent background, font size matches axis labels
axis_fs = ax.xaxis.label.get_size()
ax.legend(title='Error rate',
          ncol=2,
          fontsize=axis_fs,
          title_fontsize=axis_fs,
          loc='upper right',
          frameon=False)

# Grids
ax.grid(True, linestyle='--', linewidth=0.5)
ax.grid(True, which='minor', linestyle=':', linewidth=0.4, alpha=0.6)

plt.tight_layout()
plt.show()

