import numpy as np
import matplotlib.pyplot as plt

# Global settings
T = 300.0                          # Temperature in Kelvin
kB_kJmol = 0.0083144621           # Boltzmann constant in kJ/mol·K
kBT = kB_kJmol * T                # in kJ/mol
max_U = 8 * kBT   

def running_average(data, window_size=5):
    return np.convolve(data, np.ones(window_size)/window_size, mode='same')

import numpy as np
import matplotlib.pyplot as plt

# Globals you already had
T = 300.0
kB_kJmol = 0.0083144621
kBT = kB_kJmol * T

def running_average(data, window_size=5):
    return np.convolve(data, np.ones(window_size)/window_size, mode='same')

def invert_bond(btype,
                r_peak_min=0.25,   # ignore tiny-r bins when finding the peak
                r_peak_max=None,   # optionally cap the search window, e.g. 0.6
                smooth_win=10):

    if btype == "bond":
        bond_data = np.loadtxt("../bond.dist.new", usecols=(0,1))
        filename = "MM_bond.dat"
    elif btype == "AA":
        bond_data = np.loadtxt("AA.dist.new", usecols=(0,1))
        filename = "AA_bond.dat"
    else:
        raise ValueError("Unknown btype")

    r  = bond_data[:, 0]
    Hr = bond_data[:, 1].copy()

    # Avoid zeros
    positive = Hr > 0
    Hr[~positive] = np.min(Hr[positive]) * 1e-12

    # ---- Jacobian (r^2) with safe division; don't let super-small r blow up
    r_safe = np.clip(r, 1e-6, None)
    Pr = Hr / (r_safe**2)
    # ---------------------------------------------------------------

    # Smooth (helps stabilize the peak)
    Pr_smooth = running_average(Pr, window_size=smooth_win)

    # Boltzmann inversion
    U_r = -kBT * np.log(Pr_smooth)

    # ----- Find peak only in a trusted r-window -----
    mask = r >= r_peak_min
    if r_peak_max is not None:
        mask &= (r <= r_peak_max)

    if not np.any(mask):
        raise ValueError("Peak search window empty. Adjust r_peak_min/r_peak_max.")

    # Choose the local maximum of Pr_smooth within the mask
    # (most probable bond length ~ well minimum of U)
    local_idx = np.argmax(Pr_smooth[mask])
    # Map back to full indices
    peak_index = np.flatnonzero(mask)[local_idx]

    # Optional: refine peak by quadratic fit on a small neighborhood
    # (helps when bins are coarse)
    i0 = max(peak_index - 2, 0)
    i1 = min(peak_index + 3, len(r))  # slice is [i0, i1)
    rr = r[i0:i1]
    UU = U_r[i0:i1]
    if len(rr) >= 3:
        # Fit U(r) ~ a r^2 + b r + c; vertex gives refined minimum
        coeff = np.polyfit(rr, UU, 2)
        if coeff[0] > 0:  # parabola opens upwards
            r_vertex = -coeff[1] / (2*coeff[0])
            # Shift by U at vertex if it's inside [rr.min, rr.max], else use discrete peak
            if rr.min() <= r_vertex <= rr.max():
                U_min = np.polyval(coeff, r_vertex)
            else:
                U_min = U_r[peak_index]
        else:
            U_min = U_r[peak_index]
    else:
        U_min = U_r[peak_index]

    # Shift potential so the well minimum ~ 0 (no hard clipping)
    U_r_shifted = U_r - U_min

    # Plot
    plt.figure(figsize=(7,5))
    plt.plot(r, U_r, label="Unshifted U(r)", alpha=0.6)
    plt.plot(r, U_r_shifted, label="Shifted to min=0", lw=2)
    plt.axvline(r[peak_index], ls="--", alpha=0.6, label="peak (masked search)")
    if r_peak_max is None:
        plt.axvspan(0, r_peak_min, color='k', alpha=0.05, label="ignored for peak")
    else:
        plt.axvspan(0, r_peak_min, color='k', alpha=0.05, label="ignored for peak")
        plt.axvspan(r_peak_max, r.max(), color='k', alpha=0.05)
    plt.xlabel("Bond length (nm)")
    plt.ylabel("Potential (kJ/mol)")
    plt.legend()
    plt.grid(True)
    plt.show()

    np.savetxt(filename, np.column_stack((r, U_r_shifted)),
               fmt='%.8f', header='Bond Distance (nm)    Potential (kJ/mol)')


'''
def invert_angle(a_type):
    if a_type == "ABA":
        angle_data = np.loadtxt("../ABA.dist.new", usecols=(0,1))
        file_name = "ABA_angle.dat"
    else:
        angle_data = np.loadtxt("../BAB.dist.new", usecols=(0,1))
        file_name = "BAB_angle.dat"
    theta = angle_data[:, 0]          # in radians
    g_theta = angle_data[:, 1]

    g_theta[g_theta <= 0] = np.min(g_theta[g_theta > 0]) * 1e-6
    g_theta_smooth = running_average(g_theta, window_size=20) 
    U_theta = -kBT * np.log(g_theta_smooth)
    U_theta = np.clip(U_theta - np.min(U_theta), 0, max_U)
    #plt.plot(theta, g_theta, label="dist coarse")
    #plt.plot(theta, g_theta_smooth, label="dist smooth")
    plt.plot(theta, U_theta, label="pot")
    plt.legend()
    plt.show()
    np.savetxt(file_name, np.column_stack((theta, U_theta)), fmt='%.8f', header='Angle (rad)    Potential (kJ/mol)')

def invert_dihedral(d_type):
    if d_type == "ABAB":
        dih_data = np.loadtxt("../ABAB.dist.new", usecols=(0,1))
        file_name = "ABAB_dihedral.dat"
    else:
        dih_data = np.loadtxt("../BABA.dist.new", usecols=(0,1))
        file_name = "BABA_dihedral.dat"
    phi = dih_data[:, 0]              # in radians, assumed
    g_phi = dih_data[:, 1]

    g_phi[g_phi <= 0] = np.min(g_phi[g_phi > 0]) * 1e-6
    g_phi_smooth = running_average(g_phi, window_size=15)
    U_phi = -kBT * np.log(g_phi_smooth)
    U_phi = np.clip(U_phi - np.min(U_phi), 0, max_U)
    #plt.plot(phi, g_phi, label="coarse")
    #plt.plot(phi, g_phi_smooth, label="smooth")
    plt.plot(phi, U_phi, label="dihedral")
    plt.legend()
    plt.show()
    np.savetxt(file_name, np.column_stack((phi, U_phi)), fmt='%.8f', header='Dihedral (rad)    Potential (kJ/mol)')
'''

def invert_angle(a_type):
    if a_type == "MMM":
        angle_data = np.loadtxt("../MMM.dist.new", usecols=(0,1))
        file_name = "MMM_angle.dat"
    elif a_type == "AAB":
        angle_data = np.loadtxt("AAB.dist.new", usecols=(0,1))
        file_name = "AAB_angle.dat"
    elif a_type == "AAA":
        angle_data = np.loadtxt("AAA.dist.new", usecols=(0,1))
        file_name = "AAA_angle.dat"

    theta = angle_data[:, 0]   # radians
    H_theta = angle_data[:, 1]

    H_theta[H_theta <= 0] = np.min(H_theta[H_theta > 0]) * 1e-12

    # ---- NEW: apply Jacobian correction for angles ----
    Ptheta = H_theta / np.clip(np.sin(theta), 1e-12, None)
    # ---------------------------------------------------

    Ptheta_smooth = running_average(Ptheta, window_size=20)
    U_theta = -kBT * np.log(Ptheta_smooth)
    U_theta = np.clip(U_theta - np.min(U_theta), 0, max_U)

    # (Keep your degree conversion/interp if you like)
    theta_deg = np.degrees(theta)
    uniform_theta = np.linspace(0, 180, num=721)
    U_theta_uniform = np.interp(uniform_theta, theta_deg, U_theta)

    plt.plot(uniform_theta, U_theta_uniform, label="angle (with Jacobian)")
    plt.xlabel("Angle (degrees)"); plt.ylabel("Potential (kJ/mol)")
    plt.legend(); plt.show()

    np.savetxt(file_name, np.column_stack((uniform_theta, U_theta_uniform)),
               fmt='%.4f %.8f', header='Angle (deg)    Potential (kJ/mol)')


def invert_dihedral(d_type):
    if d_type == "MMMM":
        dih_data = np.loadtxt("../MMMM.dist.new", usecols=(0,1))
        file_name = "MMMM_dihedral.dat"
    elif d_type == "BAAA":
        dih_data = np.loadtxt("BAAA.dist.new", usecols=(0,1))
        file_name = "BAAA_dihedral.dat"
    elif d_type == "AAAB":
        dih_data = np.loadtxt("AAAB.dist.new", usecols=(0,1))
        file_name = "AAAB_dihedral.dat"
    elif d_type == "AAAA":
        dih_data = np.loadtxt("AAAA.dist.new", usecols=(0,1))
        file_name = "AAAA_dihedral.dat"
    phi_rad = dih_data[:, 0]
    g_phi = dih_data[:, 1]

    # Convert to degrees
    phi_deg = np.degrees(phi_rad)

    g_phi[g_phi <= 0] = np.min(g_phi[g_phi > 0]) * 1e-6
    g_phi_smooth = running_average(g_phi, window_size=20)

    U_phi = -kBT * np.log(g_phi_smooth)
    U_phi = np.clip(U_phi - np.min(U_phi), 0, max_U)

    # Interpolate onto uniform grid from -180 to 180 degrees
    uniform_phi = np.linspace(-180, 180, num=1441)  # 1-degree spacing
    U_phi_uniform = np.interp(uniform_phi, phi_deg, U_phi)

    plt.plot(uniform_phi, U_phi_uniform, label="Dihedral potential")
    plt.xlabel("Dihedral (degrees)")
    plt.ylabel("Potential (kJ/mol)")
    plt.legend()
    plt.show()

    np.savetxt(file_name, np.column_stack((uniform_phi, U_phi_uniform)),
               fmt='%.4f %.8f', header='Dihedral (deg)    Potential (kJ/mol)')

def invert_nonbonded(nbtype, q1, q2, rcut):
    epsilon_0 = 8.85 * 10**(-12)  # Electrostatic constant (N·m^2/C^2)
    qa = 1.602e-19 * q1  # Elementary charge (C)
    epsilon_r = 71.9  # Relative permittivity
    NA = 6.022e23  # Avogadro's number (mol^-1)
    qb = 1.602e-19 * q2
    C = (qa * qb) / (4 * 3.14159 * epsilon_0 * epsilon_r * rcut * 10 ** (-9))
    C_in_KJmol = NA * (C) / 1000
    def coulomb(x):
        V_r = NA * qa * qb * 10**(-3)/ (4 * np.pi * epsilon_r * epsilon_0 * x * 10**(-9))
        return V_r

    if nbtype=="M-NA":
        data = np.loadtxt("../M-NA.dist.new",usecols=(0,1))
        filenme="M_NA_nb.dat"
    if nbtype=="NA-NA":
        data = np.loadtxt("../NA-NA.dist.new",usecols=(0,1))
        filenme="NA_NA_nb.dat"
    if nbtype=="B-B":
        data = np.loadtxt("B-B.dist.new",usecols=(0,1))
        filenme="B_B_nb.dat"
    if nbtype=="A-NA":
        data = np.loadtxt("A-NA.dist_copy.new",usecols=(0,1))
        filenme="A_NA_nb.dat"
    if nbtype=="M-M":
        data = np.loadtxt("/home/somesh/multivalent/2d-gel/PSS-protocol/bonded_potentials/inter_range_0_6/vacuum/get_dist_1site/M-M.dist.new",usecols=(0,1))
        filenme="M_M_nb_vac.dat"

    dist = data[:,0]
    pmf = data[:,1]

    pmf[pmf <= 0] = np.min(pmf[pmf >0]) * 1e-6

    U_pmf = -kBT * np.log(pmf)
    U_pmf = np.clip(U_pmf - np.min(U_pmf), 0, max_U)
    Int_const = C_in_KJmol - U_pmf[-1]
    U_pmf[:] = U_pmf[:] + Int_const

    PMF_M = U_pmf[:] - coulomb(dist)
    np.savetxt(filenme, np.column_stack((dist, U_pmf, PMF_M)), header='r (nm)    Potential (kJ/mol)')
    plt.plot(dist, U_pmf, label=f"{nbtype}")
    plt.plot(dist, PMF_M, label="modified")
    plt.plot(dist, coulomb(dist), label="coulomb")
    plt.xlabel("r (nm)")
    plt.ylabel("Potential (KJ/mol)")
    plt.legend()
    plt.xlim([0.15, 1.3])
    plt.ylim([-10,30])
    plt.show()


invert_bond("bond")
invert_angle("MMM")
invert_dihedral("MMMM")
#invert_nonbonded("M-M", -1, -1, 0.97)
