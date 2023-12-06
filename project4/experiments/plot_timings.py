import matplotlib.pyplot as plt
import numpy as np

# Read data from file, skipping the first line
with open('data.txt', 'r') as file:
    lines = file.readlines()[1:]

# Parse the data
processors, times = zip(*[map(float, line.split()) for line in lines])

# Convert to NumPy arrays
processors = np.array(processors)
times = np.array(times)

# Create regular y-axis plot
plt.figure(figsize=(10, 5))
plt.subplot(1, 2, 1)
plt.plot(processors, times, marker='o', linestyle='-', color='darkgreen')
plt.title('Linear y axis')
plt.xlabel('Number of Processors')
plt.ylabel('Time (s)')

# Create logarithmic y-axis plot
plt.subplot(1, 2, 2)
plt.plot(processors, times, marker='o', linestyle='-', color='darkgreen')
plt.title('Logarithmic y axis')
plt.xlabel('Number of Processors')
plt.ylabel('Time (s)')
plt.yscale('log')  # Set y-axis to logarithmic scale

# Adjust boldness of frame and ticks
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['xtick.major.width'] = 1.5
plt.rcParams['ytick.major.width'] = 1.5

# Save the plot as PNG
plt.savefig('timings_plot.png')

# Display the plot
#plt.show()

