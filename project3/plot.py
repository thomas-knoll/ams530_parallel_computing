import matplotlib.pyplot as plt
import numpy as np

# Data for the speedup curves
matrix_sizes = [2**8, 2**10, 2**12]
processors = [4, 16, 64]

# Time data in seconds for each matrix size and processor count
time_data = {
    2**8: {4: 0.005759, 16: 0.004074, 64: 0.020514},
    2**10: {4: 0.315089, 16: 0.09060, 64: 0.043811},
    2**12: {4: 21.455750, 16: 5.826464, 64: 1.562086}
}

# Plotting
fig, axs = plt.subplots(len(matrix_sizes), 1, figsize=(8, 12), sharex=True)

for i, size in enumerate(matrix_sizes):
    axs[i].plot(processors, [time_data[size][p] for p in processors], marker='o', linestyle='-', label=f'Matrix Size {size}x{size}')
    axs[i].set_title(f'Speedup Curves - Matrix Size {size}x{size}')
    axs[i].set_ylabel('Time (Seconds)')
    axs[i].tick_params(axis='both', which='both', width=2, labelsize=10)

# Set common x-axis label
axs[-1].set_xlabel('Number of Processors', fontweight='bold')

plt.tight_layout(pad=3.0)
plt.subplots_adjust(top=0.94)  # Adjust the top to make space for the overall title
fig.suptitle('Speedup Curves', fontsize=16, fontweight='bold')
plt.savefig('speedup.png')
#plt.show()

