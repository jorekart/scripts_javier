import numpy as np
import matplotlib.pyplot as plt

# Read data from file
filename = 'powers.dat'  # replace with your actual file name
data = np.loadtxt(filename)

# Extract x and y columns
x = data[:, 0]
y = data[:, 1]

# Calculate dy/dx
dydx = np.gradient(y, x)

# Plot dy/dx as a function of x
plt.figure(figsize=(10, 6))
plt.plot(x, dydx, label='dy/dx', color='blue')
plt.xlabel('x')
plt.ylabel('dy/dx')
plt.title('dy/dx as a function of x')
plt.legend()
plt.grid(True)
plt.show()
