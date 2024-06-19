import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from tkinter import ttk

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from tkinter import ttk

# Define the function
def psi_n_function(psi_n, particlesource, particlesource_psin, particlesource_sig, edgeparticlesource, edgeparticlesource_psin, edgeparticlesource_sig):
    term1 = particlesource * (0.5 - 0.5 * np.tanh((psi_n - particlesource_psin) / particlesource_sig))
    term2 = edgeparticlesource * (0.5 + 0.5 * np.tanh((psi_n - edgeparticlesource_psin) / edgeparticlesource_sig))
    return term1 + term2

# Define the function to update the plot and labels
def update_plot(*args):
    particlesource = particlesource_var.get()
    particlesource_psin = particlesource_psin_var.get()
    particlesource_sig = particlesource_sig_var.get()
    edgeparticlesource = edgeparticlesource_var.get()
    edgeparticlesource_psin = edgeparticlesource_psin_var.get()
    edgeparticlesource_sig = edgeparticlesource_sig_var.get()

    psi_n = np.linspace(0, 1, 500)
    psi_n_values = psi_n_function(psi_n, particlesource, particlesource_psin, particlesource_sig, edgeparticlesource, edgeparticlesource_psin, edgeparticlesource_sig)
    
    ax.clear()
    ax.plot(psi_n, psi_n_values, label='psi_n function')
    ax.set_xlabel('psi_n')
    ax.set_ylabel('Function Value')
    ax.set_title('Plot of the psi_n function with varying parameters')
    ax.legend()
    ax.grid(True)
    canvas.draw()
    
    # Update labels
    particlesource_label.config(text=f"particlesource: {particlesource:.2f}")
    particlesource_psin_label.config(text=f"particlesource_psin: {particlesource_psin:.2f}")
    particlesource_sig_label.config(text=f"particlesource_sig: {particlesource_sig:.2f}")
    edgeparticlesource_label.config(text=f"edgeparticlesource: {edgeparticlesource:.2f}")
    edgeparticlesource_psin_label.config(text=f"edgeparticlesource_psin: {edgeparticlesource_psin:.2f}")
    edgeparticlesource_sig_label.config(text=f"edgeparticlesource_sig: {edgeparticlesource_sig:.2f}")

# Create the main window
root = tk.Tk()
root.title("psi_n Function Plotter")

# Create the figure and axis
fig, ax = plt.subplots(figsize=(10, 6))
canvas = FigureCanvasTkAgg(fig, master=root)
canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

# Create sliders
particlesource_var = tk.DoubleVar(value=1)
particlesource_psin_var = tk.DoubleVar(value=0.5)
particlesource_sig_var = tk.DoubleVar(value=0.1)
edgeparticlesource_var = tk.DoubleVar(value=1)
edgeparticlesource_psin_var = tk.DoubleVar(value=0.5)
edgeparticlesource_sig_var = tk.DoubleVar(value=0.1)

ttk.Scale(root, from_=0, to=10, orient='horizontal', variable=particlesource_var, command=update_plot).pack(fill='x', padx=10, pady=5)
particlesource_label = ttk.Label(root, text=f"particlesource: {particlesource_var.get():.2f}")
particlesource_label.pack()

ttk.Scale(root, from_=0, to=1, orient='horizontal', variable=particlesource_psin_var, command=update_plot).pack(fill='x', padx=10, pady=5)
particlesource_psin_label = ttk.Label(root, text=f"particlesource_psin: {particlesource_psin_var.get():.2f}")
particlesource_psin_label.pack()

ttk.Scale(root, from_=0.01, to=1, orient='horizontal', variable=particlesource_sig_var, command=update_plot).pack(fill='x', padx=10, pady=5)
particlesource_sig_label = ttk.Label(root, text=f"particlesource_sig: {particlesource_sig_var.get():.2f}")
particlesource_sig_label.pack()

ttk.Scale(root, from_=0, to=10, orient='horizontal', variable=edgeparticlesource_var, command=update_plot).pack(fill='x', padx=10, pady=5)
edgeparticlesource_label = ttk.Label(root, text=f"edgeparticlesource: {edgeparticlesource_var.get():.2f}")
edgeparticlesource_label.pack()

ttk.Scale(root, from_=-1, to=1, orient='horizontal', variable=edgeparticlesource_psin_var, command=update_plot).pack(fill='x', padx=10, pady=5)
edgeparticlesource_psin_label = ttk.Label(root, text=f"edgeparticlesource_psin: {edgeparticlesource_psin_var.get():.2f}")
edgeparticlesource_psin_label.pack()

ttk.Scale(root, from_=0.01, to=1, orient='horizontal', variable=edgeparticlesource_sig_var, command=update_plot).pack(fill='x', padx=10, pady=5)
edgeparticlesource_sig_label = ttk.Label(root, text=f"edgeparticlesource_sig: {edgeparticlesource_sig_var.get():.2f}")
edgeparticlesource_sig_label.pack()

# Initial plot
update_plot()

# Start the GUI event loop
root.mainloop()

