import matplotlib.pyplot as plt
import numpy as np
plt.style.use('dark_background')

# File path to the data file
data_file = 'data.o'

# Read the data
with open(data_file, 'r') as file:
    # Skip the header line
    
    # Load the rest of the data
    data = np.loadtxt(file,skiprows=8)

# Extract columns
# Assuming the columns are: Time, Pressure, Temperature
time = data[:, 0]        # First column: Time (s)
pressure = data[:, 1]    # Second column: Pressure (Pa)
temperature = data[:, 2] # Third column: Temperature (K)
parameter = data[:,1]

plt.figure(figsize=(10, 5))
plt.plot(time, parameter)
plt.title('ISP vs Time')
plt.xlabel('Time (s)')
plt.ylabel('ISP')
plt.tight_layout()
plt.show()

# # Create the first plot: Pressure vs Time
# plt.figure(figsize=(10, 5))
# plt.plot(time, pressure, label='Pressure', color='blue')
# plt.title('Pressure vs Time')
# plt.xlabel('Time (s)')
# plt.ylabel('Pressure (Pa)')
# plt.legend()
# plt.tight_layout()
# plt.savefig('pressure_vs_time.png')  # Save plot as an image
# plt.show()

# # Create the second plot: Temperature vs Time
# plt.figure(figsize=(10, 5))
# plt.plot(time, temperature, label='Temperature', color='red')
# plt.title('Temperature vs Time')
# plt.xlabel('Time (s)')
# plt.ylabel('Temperature (K)')
# plt.legend()
# plt.tight_layout()
# plt.savefig('temperature_vs_time.png')  # Save plot as an image
# plt.show()

# import matplotlib.pyplot as plt
# import numpy as np
# plt.style.use('dark_background')

# # List of data files with corresponding throat radii
# data_files = [
#     ('rt=0.0025.o', 0.0025),
#     ('rt=0.005.o', 0.005),
#     ('rt=0.0075.o', 0.0075),
#     ('rt=0.01.o', 0.01),
#     ('rt=0.0125.o', 0.0125),
#     ('rt=0.015.o', 0.015),
#     ('rt=0.0175.o', 0.0175),
#     ('rt=0.02.o', 0.02),
# ]

# # Initialize lists for legend labels
# pressure_labels = []
# temperature_labels = []

# # Initialize plots
# plt.figure(figsize=(10, 5))
# for data_file, radius in data_files:
#     with open(data_file, 'r') as file:
#         data = np.loadtxt(data_file, skiprows=8)

#     # Extract columns
#     time = data[:, 0]        # First column: Time (s)
#     pressure = data[:, 1]    # Second column: Pressure (Pa)
#     temperature = data[:, 2] # Third column: Temperature (K)

#     # Plot Pressure vs Time
#     plt.plot(time, pressure, label=f'Throat radius = {radius:.4f} m')

# plt.title('Pressure vs Time')
# plt.xlabel('Time (s)')
# plt.ylabel('Pressure (Pa)')
# plt.xlim(0, 0.5)
# plt.legend()
# plt.tight_layout()
# plt.savefig('Pressure_vs_throat_radius.png')
# plt.show()

# # Initialize plots
# plt.figure(figsize=(10, 5))
# for data_file, radius in data_files:
#     with open(data_file, 'r') as file:
#         data = np.loadtxt(data_file, skiprows=8)  # Load the rest of the data

#     # Extract columns
#     time = data[:, 0]        # First column: Time (s)
#     temperature = data[:, 2] # Third column: Temperature (K)

#     # Plot Temperature vs Time
#     plt.plot(time, temperature, label=f'Throat radius = {radius:.4f} m')

# plt.title('Temperature vs Time')
# plt.xlabel('Time (s)')
# plt.ylabel('Temperature (K)')
# plt.xlim(0, 0.5)
# plt.legend()
# plt.tight_layout()
# plt.savefig('temperature_vs_throat_radius.png')
# plt.show()