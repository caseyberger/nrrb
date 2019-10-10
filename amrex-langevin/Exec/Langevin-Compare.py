import yt
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("plotfile", type=str, help="Name of AMReX plotfile to process.")
parser.add_argument("-o", "--outfile", type=str, help="Name of output text file to write fields to.")
parser.add_argument("-c", "--compare", type=str, help="Path to a lattice save file from the non-AMReX version of NRRB.")
parser.add_argument("-csf", "--compare_sig_figs", type=int, help="Number of significant figures to compare. Default is all.")
args = parser.parse_args()

ds = yt.load(args.plotfile)
cg = ds.covering_grid(left_edge=ds.domain_left_edge, dims=ds.domain_dimensions, level=0)
coords = cg.fcoords

xs = coords[:,0]
ys = coords[:,1]
ts = coords[:,2]

dx = np.max(xs[1:] - xs[:-1])
dy = np.max(ys[1:] - ys[:-1])
dt = np.max(ts[1:] - ts[:-1])

phi1 = cg["phi1"][:,:,:].d
phi2 = cg["phi2"][:,:,:].d
phi3 = cg["phi3"][:,:,:].d
phi4 = cg["phi4"][:,:,:].d

def get_fields(i,j,k):
    return phi1[i,j,k], phi2[i,j,k], phi3[i,j,k], phi4[i,j,k]

if args.outfile:
    output_file = args.outfile
else:
    output_file = "{}_cell_values.txt".format(args.plotfile)

fout = open(output_file, "w")

fout.write("    coords   phi_1^{R}   phi_1^{I}   phi_2^{R}   phi_2^{I}\n")

for i in range(ds.domain_dimensions[0]):
    for j in range(ds.domain_dimensions[1]):
        for k in range(ds.domain_dimensions[2]):
            phi_1_Re, phi_1_Im, phi_2_Re, phi_2_Im = get_fields(i,j,k)
            fout.write("   {" + "{},{},{}".format(i,j,k) + "}" + "   {}   {}   {}   {}\n".format(phi_1_Re, phi_1_Im, phi_2_Re, phi_2_Im))
fout.close()

if args.compare:
    # Read lattice save file
    f_compare = open(args.compare, "r")

    # Read 2 header lines
    f_compare.readline()
    f_compare.readline()

    # A dictionary for holding the data
    data_compare = {}

    # Read cell & field values from file
    for l in f_compare:
        ls = l.strip().split()
        coords = ls[0][1:-1].split(',')
        x = int(coords[0])
        y = int(coords[1])
        t = int(coords[2])
        p1_Re = float(ls[1])
        p1_Im = float(ls[2])
        p2_Re = float(ls[3])
        p2_Im = float(ls[4])
        data_compare[(x,y,t)] = (p1_Re, p1_Im, p2_Re, p2_Im)

    f_compare.close()

    # Memory for norms
    InfNorm = np.zeros(4)
    L2Norm = np.zeros(4)
    InfNorm_Rel = np.zeros(4)
    L2Norm_Rel = np.zeros(4)

    # Compare the field values cell by cell
    for i in range(ds.domain_dimensions[0]):
        for j in range(ds.domain_dimensions[1]):
            for k in range(ds.domain_dimensions[2]):
                phi_1_Re_A, phi_1_Im_A, phi_2_Re_A, phi_2_Im_A = get_fields(i,j,k)
                phi_1_Re_C, phi_1_Im_C, phi_2_Re_C, phi_2_Im_C = data_compare[(i,j,k)]

                # Round the fields from the AMReX plotfile to the desired number of S.F.
                # (The previous version of NRRB uses 6 S.F. for lattice_save)
                if args.compare_sig_figs:
                    phi_1_Re_A = float('%s' % float(('%.'+'{}g').format(args.compare_sig_figs) % phi_1_Re_A))
                    phi_1_Im_A = float('%s' % float(('%.'+'{}g').format(args.compare_sig_figs) % phi_1_Im_A))
                    phi_2_Re_A = float('%s' % float(('%.'+'{}g').format(args.compare_sig_figs) % phi_2_Re_A))
                    phi_2_Im_A = float('%s' % float(('%.'+'{}g').format(args.compare_sig_figs) % phi_2_Im_A))

                # Absolute Infinity Norm
                InfNorm[0] = max(InfNorm[0], abs(phi_1_Re_A - phi_1_Re_C))
                InfNorm[1] = max(InfNorm[1], abs(phi_1_Im_A - phi_1_Im_C))
                InfNorm[2] = max(InfNorm[2], abs(phi_2_Re_A - phi_2_Re_C))
                InfNorm[3] = max(InfNorm[3], abs(phi_2_Im_A - phi_2_Im_C))

                # Absolute L2 Norm
                L2Norm[0] = L2Norm[0] + (phi_1_Re_A - phi_1_Re_C)**2
                L2Norm[1] = L2Norm[1] + (phi_1_Im_A - phi_1_Im_C)**2
                L2Norm[2] = L2Norm[2] + (phi_2_Re_A - phi_2_Re_C)**2
                L2Norm[3] = L2Norm[3] + (phi_2_Im_A - phi_2_Im_C)**2

                # Relative Infinity Norm
                InfNorm_Rel[0] = max(InfNorm_Rel[0], abs((phi_1_Re_A - phi_1_Re_C)/phi_1_Re_C))
                InfNorm_Rel[1] = max(InfNorm_Rel[1], abs((phi_1_Im_A - phi_1_Im_C)/phi_1_Im_C))
                InfNorm_Rel[2] = max(InfNorm_Rel[2], abs((phi_2_Re_A - phi_2_Re_C)/phi_2_Re_C))
                InfNorm_Rel[3] = max(InfNorm_Rel[3], abs((phi_2_Im_A - phi_2_Im_C)/phi_2_Im_C))

                # Relative L2 Norm
                L2Norm_Rel[0] = L2Norm_Rel[0] + ((phi_1_Re_A - phi_1_Re_C)/phi_1_Re_C)**2
                L2Norm_Rel[1] = L2Norm_Rel[1] + ((phi_1_Im_A - phi_1_Im_C)/phi_1_Im_C)**2
                L2Norm_Rel[2] = L2Norm_Rel[2] + ((phi_2_Re_A - phi_2_Re_C)/phi_2_Re_C)**2
                L2Norm_Rel[3] = L2Norm_Rel[3] + ((phi_2_Im_A - phi_2_Im_C)/phi_2_Im_C)**2

    # Normalize L2 Norms
    for i in range(4):
        L2Norm[i] = np.sqrt(L2Norm[i])
        L2Norm_Rel[i] = np.sqrt(L2Norm_Rel[i])

    print("Comparing {} to {} ...".format(args.plotfile, args.compare))

    if args.compare_sig_figs:
        print(" (Field values in AMReX plotfile rounded to {} significant figures for the comparison)".format(args.compare_sig_figs))

    print("  field     Abs Inf Norm   Rel Inf Norm   Abs L2 Norm    Rel L2 Norm")
    print(" phi_1_Re   {:03e}   {:03e}   {:03e}   {:03e}".format(InfNorm[0], InfNorm_Rel[0], L2Norm[0], L2Norm_Rel[0]))
    print(" phi_1_Im   {:03e}   {:03e}   {:03e}   {:03e}".format(InfNorm[1], InfNorm_Rel[1], L2Norm[1], L2Norm_Rel[1]))
    print(" phi_2_Re   {:03e}   {:03e}   {:03e}   {:03e}".format(InfNorm[2], InfNorm_Rel[2], L2Norm[2], L2Norm_Rel[2]))
    print(" phi_2_Im   {:03e}   {:03e}   {:03e}   {:03e}".format(InfNorm[3], InfNorm_Rel[3], L2Norm[3], L2Norm_Rel[3]))
