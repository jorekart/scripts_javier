"""
Template: Simple Matplotlib Plot
Author: Your Name
Description: Quick-start template for a clean, labeled plot with matplotlib.
"""
import numpy as np
import matplotlib.pyplot as plt

# -------------------
# 1. Example Data
# -------------------
x = [0, 1, 2, 3, 4, 5]
y1 = [0, 1, 4, 9, 16, 25]       # Example: y = x^2
y2 = [0, 1, 8, 27, 64, 125]     # Example: y = x^3

# -------------------
# 2. Create Figure & Axes
# -------------------
plt.figure(figsize=(8, 5))  # Width, height in inches

# -------------------
# 3. Plot Data
# -------------------
plt.plot(x, y1, label="y = x²", marker='o')
plt.plot(x, y2, label="y = x³", marker='s')

# -------------------
# 4. Customize Plot
# -------------------
plt.title("Example Plot")                 # Title
plt.xlabel("X-axis Label")                 # X-axis label
plt.ylabel("Y-axis Label")                 # Y-axis label
plt.grid(True, linestyle='--', alpha=0.6)  # Grid lines
plt.legend()                               # Show legend

# -------------------
# 5. Show or Save
# -------------------
# plt.savefig("plot.png", dpi=300)  # Save as PNG
plt.show()  # Display plot
