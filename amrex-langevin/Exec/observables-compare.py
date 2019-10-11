import yt
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-a", type=str, required=True, help="Name of AMReX observable file.")
parser.add_argument("-c", type=str, required=True, help="Name of original NRRB observable file.")
args = parser.parse_args()

# Read AMReX observable file
data_amrex = np.genfromtxt(args.a, unpack=True)

# Read original observable file
data_original = np.genfromtxt(args.c, unpack=True, skip_header=3, skip_footer=1)

# Columns
cols = ["step",
        "Re[phi^{*}phi]",
        "Im[phi^{*}phi]",
        "Re[<n>]",
        "Im[<n>]",
        "Re[<Lz>]",
        "Im[<Lz>]",
        "Re[<S>]",
        "Im[<S>]",
        "dt (sec)"]

print("Comparing {} to {} ...".format(args.a, args.c))
print("  Observable       Abs Inf Norm   Rel Inf Norm   Abs L2 Norm    Rel L2 Norm")

# Compare the observables
for i in range(1,len(cols)-2):
    da = data_amrex[i]
    do = data_original[i]
    delta_abs = da-do
    delta_rel = (da-do)/do
    InfNorm = np.linalg.norm(delta_abs, ord=np.inf)
    L2Norm = np.linalg.norm(delta_abs, ord=2)
    InfNorm_Rel = np.linalg.norm(delta_rel, ord=np.inf)
    L2Norm_Rel = np.linalg.norm(delta_rel, ord=2)
    print(" {:<15}   {:03e}   {:03e}   {:03e}   {:03e}".format(cols[i], InfNorm, InfNorm_Rel, L2Norm, L2Norm_Rel))