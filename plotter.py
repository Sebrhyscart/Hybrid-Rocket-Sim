
import matplotlib.pyplot as plt
import numpy as np

# File path to the data file
data_file = 'data.o'

# Read the data
with open(data_file, 'r') as file:
    # Skip the header line
    header = file.readline()
    
    # Load the rest of the data
    data = np.loadtxt(file)

# Extract columns
# Assuming the columns are: Time, Pressure, Temperature
time = data[:, 0]        # First column: Time (s)
pressure = data[:, 1]    # Second column: Pressure (Pa)
temperature = data[:, 2] # Third column: Temperature (K)

# Create the first plot: Pressure vs Time
plt.figure(figsize=(10, 5))
plt.plot(time, pressure, label='Pressure', color='blue')
plt.title('Pressure vs Time')
plt.xlabel('Time (s)')
plt.ylabel('Pressure (Pa)')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig('pressure_vs_time.png')  # Save plot as an image
plt.show()

# Create the second plot: Temperature vs Time
plt.figure(figsize=(10, 5))
plt.plot(time, temperature, label='Temperature', color='red')
plt.title('Temperature vs Time')
plt.xlabel('Time (s)')
plt.ylabel('Temperature (K)')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig('temperature_vs_time.png')  # Save plot as an image
plt.show()