import numpy as np

Kb = 1.38064852e-23 #* units.J / units.K  # in Joules per Kelvin
temperature = 300 #* units.K           # temperature in Kelvin (room temperature ~298.15 K)
NA = 6.02214076e23

    
def error_running(x):
    x = np.array(x)
    N = len(x)
    xmean = np.mean(x)
    xvar = np.var(x)
    acfs = []
    tau_ints = []
    k = 0
    tau_int = 0.5
    while k < 6*tau_int:
        acf = ((x[:N-k]*x[k:]).mean() - xmean*xmean) / xvar
        tau_int += acf
        N_eff = N/(2*tau_int)
        err_tau = tau_int*np.sqrt(12./N_eff)

        acfs.append(acf)
        tau_ints.append(tau_int)
        k += 1

    err_x = np.sqrt(xvar/N*2.0*tau_int)
    return xmean, xvar, err_x, tau_int, err_tau, N_eff, np.array(acfs), np.array(tau_ints)
    
def error_binning(O):
    """Block-averaging error estimate using logarithmically spaced block sizes.

    For large N this is much faster than scanning every block size from 1 to N//10,
    while still capturing the plateau in tau_int.
    """
    O = np.asarray(O, dtype=float)
    N = len(O)
    O_mean = np.mean(O)
    var_O = np.var(O, ddof=1)

    if var_O == 0.0:
        return O_mean, 0.5, N, 0.0

    # logarithmically spaced block sizes from 1 to N//10
    max_k = max(N // 10, 1)
    block_sizes = np.unique(np.geomspace(1, max_k, num=200).astype(int))

    tau_max = 0.5
    for k in block_sizes:
        N_B = N // k
        if N_B < 2:
            continue
        # vectorized: reshape and take mean of each block
        O_blocks = O[:N_B * k].reshape(N_B, k).mean(axis=1)
        var_B = np.var(O_blocks, ddof=1)
        tau_int = 0.5 * k * (var_B / var_O)
        if tau_int > tau_max:
            tau_max = tau_int

    N_eff = N / (2 * tau_max)
    var_mean = var_O / N_eff

    return O_mean, tau_max, N_eff, np.sqrt(var_mean)
