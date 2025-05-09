import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from tkinter import ttk

# Function to update the plot based on the phase difference
def update_plot(event=None):
    try:
        # Get the phase difference from the slider or entry box
        if event is None or event.widget == phase_slider:
            phase_degrees = phase_slider.get()
        else:
            phase_degrees = float(phase_entry.get())
            phase_slider.set(phase_degrees)  # Sync the slider with the entry box

        # Convert degrees to radians for calculations
        phase_radians = np.radians(phase_degrees)

        # Clear the previous plot
        ax.clear()

        # Parameters
        Ex_magnitude = 1.0  # Amplitude of Ex
        Ey_magnitude = 1.0  # Amplitude of Ey
        time = np.linspace(0, 2 * np.pi, 1000)  # Time array

        # Electric field components
        Ex = Ex_magnitude * np.cos(time)  # Ex = |Ex| * cos(ωt)
        Ey = Ey_magnitude * np.cos(time + phase_radians)  # Ey = |Ey| * cos(ωt + Δϕ)

        # Plotting the polarization ellipse
        ax.plot(Ex, Ey, label=f"Phase Difference = {phase_degrees:.1f}°")
        ax.set_title("Polarization State of Light")
        ax.set_xlabel("Ex (Horizontal)")
        ax.set_ylabel("Ey (Vertical)")
        ax.grid(True)
        ax.axhline(0, color='black', linewidth=0.5)
        ax.axvline(0, color='black', linewidth=0.5)
        ax.legend()
        ax.set_xlim(-1.5, 1.5)
        ax.set_ylim(-1.5, 1.5)

        # Add text to indicate the rotation direction
        if phase_degrees % 360 != 0 and phase_degrees % 360 != 180:  # Not linear polarization
            if phase_degrees % 360 < 180:  # Right-handed rotation
                ax.text(0.5, 1.1, "Right-Handed Rotation", color='red', fontsize=12, ha='center',
                        transform=ax.transAxes)
            else:  # Left-handed rotation
                ax.text(0.5, 1.1, "Left-Handed Rotation", color='blue', fontsize=12, ha='center',
                        transform=ax.transAxes)
        else:
            ax.text(0.5, 1.1, "No Rotation (Linear Polarization)", color='black', fontsize=12, ha='center',
                    transform=ax.transAxes)

        # Redraw the canvas
        canvas.draw()
    except ValueError:
        tk.messagebox.showerror("Error", "Please enter a valid number for the phase difference.")

# Create the main window
root = tk.Tk()
root.title("Polarization Visualization")
plot_frame = ttk.Frame(root)
plot_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

# Create a matplotlib figure and axis
fig, ax = plt.subplots(figsize=(6, 6))
canvas = FigureCanvasTkAgg(fig, master=plot_frame)
canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
input_frame = ttk.Frame(root)
input_frame.pack(side=tk.BOTTOM, fill=tk.X)

# Add a slider for the phase difference (in degrees)
ttk.Label(input_frame, text="Phase Difference (°):").pack(side=tk.LEFT, padx=5, pady=5)
phase_slider = tk.Scale(input_frame, from_=0, to=360, resolution=0.1, orient=tk.HORIZONTAL, length=300)
phase_slider.set(0)
phase_slider.pack(side=tk.LEFT, padx=5, pady=5)

# For entry
phase_entry = ttk.Entry(input_frame, width=10)
phase_entry.pack(side=tk.LEFT, padx=5, pady=5)
phase_entry.insert(0, "0")

# Bind the slider and entry to the update_plot function
phase_slider.bind("<Motion>", update_plot)
phase_entry.bind("<Return>", update_plot)

update_plot()
root.mainloop()