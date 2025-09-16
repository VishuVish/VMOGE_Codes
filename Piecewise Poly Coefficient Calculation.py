import numpy as np

def labview_style_coeffs(coeffs_desc):
    """Convert NumPy-style coefficients (highest degree first) to LabVIEW-style (ascending)."""
    return coeffs_desc[::-1]

def generate_piecewise_polynomials_auto_axis(file_name, low_thresh_mT=80, degree_low=3, degree_high=1):
    """
    Generate piecewise inverse polynomials for LabVIEW with automatic swept-axis detection.
    Units: everything converted to Tesla before fitting.
    Skips missing/NaN values in the data.
    
    Parameters:
        file_name : str
            Text file with columns: Applied (mT), Bx (mT), By (mT), Bz (mT)
        low_thresh_mT : float
            Threshold to separate low/high field regions (in mT, will be converted to Tesla internally)
        degree_low : int
            Polynomial degree for low-field region
        degree_high : int
            Polynomial degree for high-field region
    """
    # Load data
    data = np.loadtxt(file_name)
    applied_mT = data[:,0]
    Bx_mT, By_mT, Bz_mT = data[:,1], data[:,2], data[:,3]

    # --- Convert all to Tesla ---
    applied = applied_mT / 1000.0
    Bx, By, Bz = Bx_mT / 1000.0, By_mT / 1000.0, Bz_mT / 1000.0
    low_thresh = low_thresh_mT / 1000.0  # threshold in Tesla

    # --- Detect swept axis ---
    stds = [np.nanstd(Bx), np.nanstd(By), np.nanstd(Bz)]  # ignore NaNs in std calculation
    axes = ['X', 'Y', 'Z']
    swept_idx = np.nanargmax(stds)
    swept_axis = axes[swept_idx]

    # Select measured values for the swept axis
    measured_dict = {'X': Bx, 'Y': By, 'Z': Bz}
    measured = measured_dict[swept_axis]

    print(f"Detected swept axis: {swept_axis}")

    # --- Remove NaNs from applied & measured ---
    valid_mask = ~np.isnan(measured) & ~np.isnan(applied)
    applied = applied[valid_mask]
    measured = measured[valid_mask]

    # --- Split into low/high-field regions ---
    low_mask = np.abs(measured) <= low_thresh
    high_mask = np.abs(measured) > low_thresh

    measured_low = measured[low_mask]
    applied_low = applied[low_mask]

    measured_high = measured[high_mask]
    applied_high = applied[high_mask]

    # --- Fit inverse polynomials if data exists ---
    if len(measured_low) > 0:
        coeffs_low = np.polyfit(measured_low, applied_low, degree_low)
        coeffs_low_labview = labview_style_coeffs(coeffs_low)
        print("\nLow-field inverse polynomial (Measured -> Applied):")
        print("Degree:", degree_low)
        print("NumPy order:", coeffs_low)
        print("LabVIEW order:", coeffs_low_labview)
    else:
        coeffs_low_labview = None
        print("\nNo valid low-field data to fit.")

    if len(measured_high) > 0:
        coeffs_high = np.polyfit(measured_high, applied_high, degree_high)
        coeffs_high_labview = labview_style_coeffs(coeffs_high)
        print("\nHigh-field inverse polynomial (Measured -> Applied):")
        print("Degree:", degree_high)
        print("NumPy order:", coeffs_high)
        print("LabVIEW order:", coeffs_high_labview)
    else:
        coeffs_high_labview = None
        print("\nNo valid high-field data to fit.")

    return swept_axis, coeffs_low_labview, coeffs_high_labview

# --- Example usage ---
file_name = "Measured_Mag_-Y.txt"  # replace with your file
swept_axis, coeffs_low, coeffs_high = generate_piecewise_polynomials_auto_axis(file_name)
