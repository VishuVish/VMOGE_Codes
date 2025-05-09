import numpy as np
import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# Constants
alpha = 7.31  # A/T, depends on the magnet design


# Conversion functions
def polar_to_cartesian(H, phi_m, theta_m):
    H_x = H * np.sin(np.radians(theta_m)) * np.cos(np.radians(phi_m))
    H_y = H * np.sin(np.radians(theta_m)) * np.sin(np.radians(phi_m))
    H_z = H * np.cos(np.radians(theta_m))
    return H_x, H_y, H_z


def cartesian_to_polar(H_x, H_y, H_z):
    H = np.sqrt(H_x ** 2 + H_y ** 2 + H_z ** 2)
    phi_m = np.degrees(np.arctan2(H_y, H_x))
    theta_m = np.degrees(np.arccos(H_z / H))
    return H, phi_m, theta_m


def currents_to_cartesian(I1, I2, I3, I4):
    H_x = (-I1 - I2 + I3 + I4) / (4 * alpha)
    H_y = (I1 + I2 + I3 + I4) / (4 * alpha)
    H_z = (I1 - I2 - I3 + I4) / (4 * alpha)
    return H_x, H_y, H_z

def normalize_magnetic_field(H_x, H_y, H_z, H):
    magnitude = np.sqrt(H_x**2 + H_y**2 + H_z**2)
    if magnitude > 0:  # Avoid division by zero
        H_x_norm = H_x / magnitude
        H_y_norm = H_y / magnitude
        H_z_norm = H_z / magnitude
    else:
        H_x_norm, H_y_norm, H_z_norm = 0, 0, 0  # Default to zero if magnitude is zero
    # Scale to desired magnitude
    H_x_final = H * H_x_norm
    H_y_final = H * H_y_norm
    H_z_final = H * H_z_norm
    return H_x_final, H_y_final, H_z_final


def cartesian_to_currents(H_x_final, H_y_final, H_z_final):
    I1 = 2 * alpha * H_y_final - 2 * alpha * H_z_final
    I2 = 2 * alpha * H_x_final + 2 * alpha * H_z_final
    I3 = -2 * alpha * H_x_final + 2 * alpha * H_y_final
    I4 = 0  # Assuming minimum norm solution
    return I1, I2, I3, I4

def determine_vmoge_mode(H_x, H_y, H_z):
    # Check for LP-VMOGE (x-z plane, H_y = 0)
    if np.isclose(H_y, 0, atol=1e-6) and not np.isclose(H_x, 0, atol=1e-6) and not np.isclose(H_z, 0, atol=1e-6):
        return "LP-VMOGE"

    # Check for TP-VMOGE (y-z plane, H_x = 0)
    elif np.isclose(H_x, 0, atol=1e-6) and not np.isclose(H_y, 0, atol=1e-6) and not np.isclose(H_z, 0, atol=1e-6):
        return "TP-VMOGE"

    # Check for LT-VMOGE (x-y plane, H_z = 0)
    elif np.isclose(H_z, 0, atol=1e-6) and not np.isclose(H_x, 0, atol=1e-6) and not np.isclose(H_y, 0, atol=1e-6):
        return "LT-VMOGE"

    # Default to 3D-VMOGE if none of the above conditions are met
    else:
        return "3D-VMOGE"

# Function to calculate light path
def calculate_light_path(gamma):
    # Incident light direction (from above the sample, pointing to (0, 0, 0))
    gamma_rad = np.radians(gamma)
    incident_dir = np.array([np.sin(gamma_rad), 0, -np.cos(gamma_rad)])  # Light comes from +z direction at angle gamma
    # Reflected light direction (assuming perfect reflection)
    reflected_dir = np.array([np.sin(gamma_rad), 0, np.cos(gamma_rad)])
    return incident_dir, reflected_dir

def update(event=None):
    try:
        if mode.get() == "Polar":
            # Fixed magnitude from H entry
            H = float(h_entry.get())
            # Orientation angles from sliders or entry fields
            phi_m = phi_var.get()
            theta_m = theta_var.get()
            # Calculate Cartesian components
            H_x, H_y, H_z = polar_to_cartesian(H, phi_m, theta_m)
            # Calculate currents
            I1, I2, I3, I4 = cartesian_to_currents(H_x, H_y, H_z)

        elif mode.get() == "Cartesian":
            # Cartesian components from entry fields
            H_x = float(cartesian_entries[0].get())
            H_y = float(cartesian_entries[1].get())
            H_z = float(cartesian_entries[2].get())

            # Recalculate magnitude and orientation angles
            H, phi_m, theta_m = cartesian_to_polar(H_x, H_y, H_z)
            # Calculate currents
            I1, I2, I3, I4 = cartesian_to_currents(H_x, H_y, H_z)

        elif mode.get() == "Currents":
            # Only update if values are manually changed
            I1 = float(current_entries[0].get())
            I2 = float(current_entries[1].get())
            I3 = float(current_entries[2].get())
            I4 = float(current_entries[3].get())

            # Recalculate Cartesian components from currents
            H_x, H_y, H_z = currents_to_cartesian(I1, I2, I3, I4)

            # Check the toggle state to decide whether to normalize
            if normalize_field.get():  # Field Orbit Mode (Normalize)
                H = float(h_entry.get())  # Fixed magnitude from polar input
                H_x, H_y, H_z = normalize_magnetic_field(H_x, H_y, H_z, H)
            else:  # Hysteresis Mode (Do not normalize)
                H, phi_m, theta_m = cartesian_to_polar(H_x, H_y, H_z)
                # Update slider positions
                phi_var.set(phi_m)
                theta_var.set(theta_m)

        # Update all entry fields
        h_entry.delete(0, tk.END)
        h_entry.insert(0, f"{H:.2f}")

        cartesian_entries[0].delete(0, tk.END)
        cartesian_entries[0].insert(0, f"{H_x:.2f}")
        cartesian_entries[1].delete(0, tk.END)
        cartesian_entries[1].insert(0, f"{H_y:.2f}")
        cartesian_entries[2].delete(0, tk.END)
        cartesian_entries[2].insert(0, f"{H_z:.2f}")

        current_entries[0].delete(0, tk.END)
        current_entries[0].insert(0, f"{I1:.2f}")
        current_entries[1].delete(0, tk.END)
        current_entries[1].insert(0, f"{I2:.2f}")
        current_entries[2].delete(0, tk.END)
        current_entries[2].insert(0, f"{I3:.2f}")
        current_entries[3].delete(0, tk.END)
        current_entries[3].insert(0, f"{I4:.2f}")

        # Determine the VMOGE mode using Cartesian coordinates
        vmoge_mode = determine_vmoge_mode(H_x, H_y, H_z)
        mode_label.config(text=f"VMOGE Mode: {vmoge_mode}")

        # Get light angle
        gamma = float(light_angle_entry.get())
        incident_dir, reflected_dir = calculate_light_path(gamma)

        # Update the plot
        ax.clear()

        # Plot the magnetic field vector
        ax.quiver(0, 0, 0, H_x, H_y, H_z, color='r', label='Applied Magnetic Field Direction')

        # Plot incident and reflected light
        start_x = -incident_dir[0]
        start_y = -incident_dir[1]
        start_z = -incident_dir[2]
        ax.quiver(start_x, start_y, start_z, -start_x, -start_y, -start_z, color='green', label='Incident Light')
        ax.quiver(0, 0, 0, reflected_dir[0], reflected_dir[1], reflected_dir[2], color='orange', label='Reflected Light')

        # Plot the sample surface (x-y plane)
        xx, yy = np.meshgrid(np.linspace(-1, 1, 2), np.linspace(-1, 1, 2))
        zz = np.zeros_like(xx)
        ax.plot_surface(xx, yy, zz, color='gray', alpha=0.5, label='Sample Surface (x-y plane)')

        # Set plot limits and labels
        ax.set_xlim([-1, 1])
        ax.set_ylim([-1, 1])
        ax.set_zlim([-1, 1])

        # Set specific tick values for x, y, and z axes
        ax.set_xticks([-1, -0.5, 0, 0.5, 1])
        ax.set_yticks([-1, -0.5, 0, 0.5, 1])
        ax.set_zticks([-1, -0.5, 0, 0.5, 1])

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.legend(loc='upper right')
        canvas.draw()

    except ValueError as e:
        print(f"Error in update function: {e}")
        pass

#Hysteresis or Field Orbit
def handle_mode_change():
    if mode.get() == "Polar":
        # Disable the toggle button in Polar Mode
        toggle_button.config(state=tk.DISABLED)
    else:
        # Enable the toggle button in Cartesian/Current Mode
        toggle_button.config(state=tk.NORMAL)
    # Update the calculations
    update()

def toggle_measurement_mode():
    if normalize_field.get():  # Toggle is on (Field Orbit Mode)
        # Switch to Cartesian or Current Mode
        if mode.get() != "Cartesian" and mode.get() != "Currents":
            mode.set("Cartesian")  # Default to Cartesian Mode
    else:  # Toggle is off (Hysteresis Mode)
        # No need to switch modes
        pass

    # Update the calculations
    update()

# UI setup
root = tk.Tk()
root.title("Magnetic Field Converter")

# Mode selection
mode = tk.StringVar(value="Polar")
mode_frame = ttk.LabelFrame(root, text="Input Mode")
mode_frame.grid(row=0, column=0, padx=10, pady=10, sticky="ew")
ttk.Radiobutton(mode_frame, text="Polar", variable=mode, value="Polar").grid(row=0, column=0,
                                                                                             sticky="w")
ttk.Radiobutton(mode_frame, text="Cartesian", variable=mode, value="Cartesian").grid(row=1, column=0,
                                                                                                     sticky="w")
ttk.Radiobutton(mode_frame, text="Currents", variable=mode, value="Currents").grid(row=2, column=0,
                                                                                                   sticky="w")

# Bind the mode change function to the mode variable
mode.trace_add("write", lambda *args: handle_mode_change())

# Add a toggle button for Hysteresis/Field Orbit mode
toggle_frame = ttk.LabelFrame(root, text="Measurement Mode")
toggle_frame.grid(row=6, column=0, padx=10, pady=10, sticky="ew")

# Boolean variable to track the toggle state
normalize_field = tk.BooleanVar(value=False)  # False = Hysteresis, True = Field Orbit

# Toggle button (Checkbutton)
toggle_button = ttk.Checkbutton(
    toggle_frame,
    text="Field Orbit Mode (Normalize Field)",
    variable=normalize_field,
    command=toggle_measurement_mode  # Call toggle_measurement_mode when the toggle state changes
)
toggle_button.grid(row=0, column=0, sticky="w")

# Polar entries
polar_frame = ttk.LabelFrame(root, text="Polar Coordinates")
polar_frame.grid(row=1, column=0, padx=10, pady=10, sticky="ew")

# H field entry
ttk.Label(polar_frame, text="H (T)").grid(row=0, column=0, sticky="w")
h_entry = ttk.Entry(polar_frame)
h_entry.grid(row=0, column=1, sticky="ew")
h_entry.insert(0, "1.0")  # Default value

# Phi_m slider
ttk.Label(polar_frame, text="φ_m (°)").grid(row=1, column=0, sticky="w")
phi_var = tk.DoubleVar(value=0)  # Default value of 0
phi_slider = ttk.Scale(polar_frame, from_=0, to=360, orient="horizontal", variable=phi_var, command=lambda _: update())
phi_slider.grid(row=1, column=1, sticky="ew")
phi_value_label = ttk.Label(polar_frame, text="0")
phi_entry = ttk.Entry(polar_frame, width=5)
phi_entry.grid(row=1, column=2, padx=5)
phi_entry.insert(0, "0")  # Default value


# Theta_m slider
ttk.Label(polar_frame, text="θ_m (°)").grid(row=2, column=0, sticky="w")
theta_var = tk.DoubleVar(value=90)  # Default value of 90
theta_slider = ttk.Scale(polar_frame, from_=0, to=360, orient="horizontal", variable=theta_var,
                         command=lambda _: update())
theta_slider.grid(row=2, column=1, sticky="ew")
theta_value_label = ttk.Label(polar_frame, text="90")
theta_entry = ttk.Entry(polar_frame, width=5)
theta_entry.grid(row=2, column=2, padx=5)
theta_entry.insert(0, "90")  # Default value



# Update entry fields when sliders change
def update_angle_labels(*args):
    phi_entry.delete(0, tk.END)
    phi_entry.insert(0, f"{int(phi_var.get())}")
    theta_entry.delete(0, tk.END)
    theta_entry.insert(0, f"{int(theta_var.get())}")

# Register callbacks
phi_var.trace_add("write", update_angle_labels)
theta_var.trace_add("write", update_angle_labels)

# Cartesian entries
cartesian_frame = ttk.LabelFrame(root, text="Cartesian Coordinates")
cartesian_frame.grid(row=2, column=0, padx=10, pady=10, sticky="ew")
cartesian_labels = ["H_x (T)", "H_y (T)", "H_z (T)"]
cartesian_entries = []
for i, label in enumerate(cartesian_labels):
    ttk.Label(cartesian_frame, text=label).grid(row=i, column=0, sticky="w")
    entry = ttk.Entry(cartesian_frame)
    entry.grid(row=i, column=1, sticky="ew")
    cartesian_entries.append(entry)

# Current entries
current_frame = ttk.LabelFrame(root, text="Currents (A)")
current_frame.grid(row=3, column=0, padx=10, pady=10, sticky="ew")
current_labels = ["I1", "I2", "I3", "I4"]
current_entries = []
for i, label in enumerate(current_labels):
    ttk.Label(current_frame, text=label).grid(row=i, column=0, sticky="w")
    entry = ttk.Entry(current_frame)
    entry.grid(row=i, column=1, sticky="ew")
    current_entries.append(entry)

# Light angle entry
light_angle_frame = ttk.LabelFrame(root, text="Light Angle (γ)")
light_angle_frame.grid(row=4, column=0, padx=10, pady=10, sticky="w")
ttk.Label(light_angle_frame, text="γ (°)").grid(row=0, column=0, sticky="w")
light_angle_entry = ttk.Entry(light_angle_frame)
light_angle_entry.grid(row=0, column=1, sticky="ew")
light_angle_entry.insert(0, "45")  # Default angle

# Bind entry fields to update only on manual input
h_entry.bind("<Return>", update)
phi_entry.bind("<Return>", update)
theta_entry.bind("<Return>", update)
light_angle_entry.bind("<Return>", update)

for entry in cartesian_entries:
    entry.bind("<Return>", update)

for entry in current_entries:
    entry.bind("<Return>", update)

# Bind entry fields to update sliders
phi_entry.bind("<Return>", lambda event: phi_var.set(float(phi_entry.get())))
theta_entry.bind("<Return>", lambda event: theta_var.set(float(theta_entry.get())))

# Update sliders when entry fields change
phi_entry.bind("<FocusOut>", lambda event: phi_var.set(float(phi_entry.get())))
theta_entry.bind("<FocusOut>", lambda event: theta_var.set(float(theta_entry.get())))

# VMOGE mode label
mode_label = ttk.Label(root, text="VMOGE Mode: Unknown")
mode_label.grid(row=0, column=2, padx=10, pady=10)

# Update button
update_button = ttk.Button(root, text="Update", command=update)
update_button.grid(row=5, column=0, padx=10, pady=10)

# Plot
fig = plt.figure(figsize=(11, 7))
ax = fig.add_subplot(111, projection='3d')
canvas = FigureCanvasTkAgg(fig, master=root)
canvas.get_tk_widget().grid(row=0, column=1, rowspan=6, padx=10, pady=10)

handle_mode_change()

# Initial update
update()
root.mainloop()