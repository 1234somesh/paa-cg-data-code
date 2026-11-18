import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def compute_force(positions, pmf_values):
    """
    Compute the force as the negative derivative of the PMF using the central finite difference method.

    Parameters:
    positions (numpy array): Array of position values.
    pmf_values (numpy array): Array of PMF values corresponding to the positions.

    Returns:
    forces (numpy array): Estimated force at each position.
    """
    forces = np.zeros_like(pmf_values)

    # Central difference for interior points
    for i in range(1, len(positions) - 1):
        forces[i] = -(pmf_values[i+1] - pmf_values[i-1]) / (positions[i+1] - positions[i-1])

    # Forward difference for the first point
    forces[0] = -(pmf_values[1] - pmf_values[0]) / (positions[1] - positions[0])

    # Backward difference for the last point
    forces[-1] = -(pmf_values[-1] - pmf_values[-2]) / (positions[-1] - positions[-2])

    return forces

def harmonic(r, k, r0, V0):
    return 0.5 * k * (r - r0)**2 + V0

def extrapolate_pmf(n_points=4):
    data = np.loadtxt('MM_bond.dat')
    positions = data[242:472,0]
    pmf_values = data[242:472,1]
    # Left-side fit
    r_left = positions[:n_points]
    V_left = pmf_values[:n_points]
    popt_left, _ = curve_fit(harmonic, r_left, V_left, p0=[1000, r_left[np.argmin(V_left)], np.min(V_left)], maxfev=1000000)

    # Right-side fit
    r_right = positions[-n_points:]
    V_right = pmf_values[-n_points:]
    popt_right, _ = curve_fit(harmonic, r_right, V_right, p0=[1000, r_right[np.argmin(V_right)], np.min(V_right)], maxfev=1000000)

    # Extrapolate beyond the current range
    dr = np.mean(np.diff(positions))
    r_min = positions[0] - 0.24
    r_max = positions[-1] + 1.50

    r_extra_left = np.arange(r_min, positions[0], dr)
    r_extra_right = np.arange(positions[-1] + dr, r_max + dr, dr)

    V_extra_left = harmonic(r_extra_left, *popt_left)
    V_extra_right = harmonic(r_extra_right, *popt_right)

    # Concatenate all
    positions_full = np.concatenate([r_extra_left, positions, r_extra_right])
    pmf_full = np.concatenate([V_extra_left, pmf_values, V_extra_right])
    forces = compute_force(positions_full, pmf_full)
    plt.plot(positions, pmf_values, label="original")
    plt.plot(positions_full, pmf_full, label="extrapolated")
    plt.legend()
    plt.show()
    return positions_full, pmf_full, forces

def write_xvg(filename, positions, pmf_values, forces):
    with open(filename, 'w') as f:
        f.write("@    title \"Tabulated Nonbonded Interaction\"\n")
        f.write("@    xaxis  label \"Position (nm)\"\n")
        f.write("@    yaxis  label \"PMF (kJ/mol)\"\n")
        f.write("@TYPE xy\n")
        for i in range(len(positions)):
            row = [positions[i]] + [pmf_values[i], forces[i]]
            f.write("  ".join(f"{val: .6f}" for val in row) + "\n")

def main():
    #extended_positions, extended_pmf, forces = extrapolate_pmf(8)
    #write_xvg("table_b0.xvg", extended_positions, extended_pmf, forces)
    #theta, potential = np.loadtxt("MMM_angle.dat", unpack=True)
    #forces = compute_force(theta, potential)
    #write_xvg("table_a0.xvg", theta, potential, forces)
    theta, potential = np.loadtxt("MMMM_dihedral.dat", unpack=True)
    forces = compute_force(theta, potential)
    write_xvg("table_d0.xvg", theta, potential, forces)
    print(f"Done: Generated")

if __name__ == "__main__":
    main()
