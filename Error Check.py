import numpy as np
import matplotlib.pyplot as plt

def analyze_file(fname, tolerance_pct=3.0):
    # Load data
    data = np.loadtxt(fname)
    Applied_file = data[:,0]
    Bx, By, Bz = data[:,1], data[:,2], data[:,3]
    n = len(Applied_file)

    # ---- Ask user if they want to override applied values ----
    override = input("Override applied field values from file? (y/n): ").strip().lower()
    if override == "y":
        start = float(input("Enter magnitude start (mT): "))
        end   = float(input("Enter magnitude end (mT): "))
        step  = float(input("Enter step size (mT): "))
        Applied = np.arange(start, end + step, step)
        if len(Applied) != n:
            print(f"⚠️ Warning: Applied length ({len(Applied)}) "
                  f"≠ data length ({n}). Truncating to match.")
            Applied = Applied[:n]
    else:
        Applied = Applied_file.copy()

    # Detect applied axis (largest variation in measured fields)
    ranges = [np.ptp(Bx), np.ptp(By), np.ptp(Bz)]
    axes = ['x', 'y', 'z']
    applied_axis_idx = np.argmax(ranges)
    applied_axis = axes[applied_axis_idx]

    # Grab the measured applied axis
    meas_dict = {'x': Bx, 'y': By, 'z': Bz}
    meas_applied = meas_dict[applied_axis]

    # --- Robust sign alignment (fix for the 200% error) ---
    def dominant_sign(arr):
        """Return a robust sign for an array: use mean if non-zero, else slope."""
        if np.isnan(arr).all():
            return 1
        m = np.nanmean(arr)
        if abs(m) > 1e-12:
            return np.sign(m)
        # fallback to slope over index (handles symmetric arrays)
        p = np.polyfit(np.arange(len(arr)), arr, 1)
        return np.sign(p[0]) if abs(p[0]) > 1e-12 else 1

    meas_sign = dominant_sign(meas_applied)
    applied_sign = dominant_sign(Applied)

    # Only flip Applied when its overall sign disagrees with the measured axis sign.
    if applied_sign != meas_sign:
        Applied = -Applied

    # Determine sweep direction from measured axis
    sweep_dir = "positive" if meas_sign > 0 else "negative"

    # Make plotting x-axis start at 0
    applied = Applied.copy()
    x_axis = applied - applied[0]

    print(f"Detected applied axis: {applied_axis.upper()} ({sweep_dir} sweep)")

    # Assign the other axes and labels
    if applied_axis == 'x':
        meas_other1, meas_other2 = By, Bz
        labels_signed = ["Error Bx (%)", "Error By (%) leakage", "Error Bz (%) leakage"]
    elif applied_axis == 'y':
        meas_other1, meas_other2 = Bx, Bz
        labels_signed = ["Error By (%)", "Error Bx (%) leakage", "Error Bz (%) leakage"]
    else:  # z
        meas_other1, meas_other2 = Bx, By
        labels_signed = ["Error Bz (%)", "Error Bx (%) leakage", "Error By (%) leakage"]

    # --- Percent errors with safe division (mask zeros) ---
    mask_nonzero = applied != 0
    if not np.any(mask_nonzero):
        raise ValueError("All applied values are zero — percent error undefined.")

    # initialize arrays with NaN for undefined (applied==0)
    err_applied = np.full(n, np.nan)
    err_leak1   = np.full(n, np.nan)
    err_leak2   = np.full(n, np.nan)

    err_applied[mask_nonzero] = (meas_applied[mask_nonzero] - applied[mask_nonzero]) / applied[mask_nonzero] * 100.0
    err_leak1[mask_nonzero]   = meas_other1[mask_nonzero] / applied[mask_nonzero] * 100.0
    err_leak2[mask_nonzero]   = meas_other2[mask_nonzero] / applied[mask_nonzero] * 100.0

    errors_signed = [err_applied, err_leak1, err_leak2]
    errors_abs = [np.abs(err) for err in errors_signed]
    labels_abs = [f"|{lab}|" for lab in labels_signed]

    # -------- Summary stats (ignore NaNs) --------
    def summarize(name, signed, absolute):
        valid = ~np.isnan(absolute)
        count_valid = np.sum(valid)
        if count_valid == 0:
            print(f"{name}: no valid (non-zero applied) points to compute percent error.")
            return
        mean_signed = np.nanmean(signed)
        mean_abs = np.nanmean(absolute)
        max_abs = np.nanmax(absolute)
        frac_within_tol = np.sum(absolute[valid] <= tolerance_pct) / count_valid
        print(f"{name}: mean_signed={mean_signed:.2f}%, mean_abs={mean_abs:.2f}%, "
              f"max_abs={max_abs:.2f}%, within ±{tolerance_pct}%={frac_within_tol*100:.1f}% "
              f"({count_valid}/{n} valid points)")

    print("\nSummary (percent):")
    for lab, s, a in zip(labels_signed, errors_signed, errors_abs):
        summarize(lab, s, a)

    if np.any(applied == 0):
        print("⚠️ Note: Some applied values are exactly zero — percent error at those points is undefined and omitted.")

    # -------- Plot 1: Signed errors --------
    plt.figure(figsize=(11,6))
    markers = ['^', 'o', 's']
    linestyles = ['-', '--', '--']
    for err, lab, m, ls in zip(errors_signed, labels_signed, markers, linestyles):
        plt.plot(x_axis, err, label=lab, marker=m, markersize=3, linestyle=ls)

    plt.fill_between(x_axis, -tolerance_pct, tolerance_pct,
                     color='gray', alpha=0.15, label=f"±{tolerance_pct}% tolerance")
    plt.axhline(0, color='k', lw=0.8)
    plt.xlabel(f"Applied Magnetic Field along {applied_axis.upper()} mT")
    plt.ylabel("Recorded Error (%)")
    plt.title(f"Recorded Percentage Error with Correction — Applied along {applied_axis.upper()} axis ({sweep_dir} sweep)")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    if sweep_dir == "negative":
        plt.gca().invert_xaxis()
    plt.show()

    # -------- Plot 2: Absolute errors --------
    plt.figure(figsize=(11,6))
    for err, lab, m, ls in zip(errors_abs, labels_abs, markers, linestyles):
        plt.plot(x_axis, err, label=lab, marker=m, markersize=3, linestyle=ls)

    plt.axhline(tolerance_pct, color='red', linestyle=':', lw=1.5,
                label=f"{tolerance_pct}% tolerance")
    plt.xlabel(f"Applied Magnetic Field along {applied_axis.upper()} mT")
    plt.ylabel("Absolute Error (%)")
    plt.title(f"Absolute Percentage Error with Correction — Applied along {applied_axis.upper()} axis ({sweep_dir} sweep)")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    if sweep_dir == "negative":
        plt.gca().invert_xaxis()
    plt.show()

if __name__ == "__main__":
    analyze_file("Measured_Mag_-Y_Final.txt", tolerance_pct=3.0)
