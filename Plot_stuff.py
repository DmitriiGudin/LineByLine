from __future__ import division
import numpy as np
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


vel_file_list = "Svel.list"
corr_file = "RAVE_rv_corr.txt"


def get_predicted_vels (File):
    vel, vel_err = [], []
    for f in np.genfromtxt(File, dtype=str, delimiter='\t'):
        with open(f, 'r') as f_i:
            lines = f_i.read().splitlines()
            vel.append(float(lines[-1].split()[-2]))
            vel_err.append(float(lines[-1].split()[-1]))
    return np.array(vel), np.array(vel_err)


def get_RAVE_vels_corr (File):
    vel, vel_n, vel_corr = [], [], []
    with open(File, 'r') as f:
        lines = f.read().splitlines()
        for l in lines:
            vel.append(float(l.split()[1]))
            vel_n.append(float(l.split()[3]))
            vel_corr.append(float(l.split()[-1]))
    return np.array(vel), np.array(vel_n), np.array(vel_corr)


def plot_stuff (vel0, vel0_n, vel1):
    plt.clf()
    plt.title("Line-by-line RV calculation", size=24)
    plt.xlabel("Radial velocity (km/s)", size=24)
    plt.ylabel("Radial velocity residual (km/s)", size=24)
    plt.tick_params(labelsize=18)
    plt.scatter(vel1, vel1-vel0, color='black', s=20)
    plt.gcf().set_size_inches(25.6, 14.4)
    plt.gcf().savefig("RV_vs_res.png", dpi=100)
    plt.close()

    plt.clf()
    plt.title("Line-by-line RV calculation", size=24)
    plt.xlabel("Radial velocity residual(km/s)", size=24)
    plt.ylabel("Number of occurrences", size=24)
    plt.tick_params(labelsize=18)
    plt.hist(vel1-vel0, linewidth=2, color='black', histtype='step', bins=50)
    plt.gcf().set_size_inches(25.6, 14.4)
    plt.gcf().savefig("RV_res.png", dpi=100)
    plt.close()


if __name__ == '__main__':
    vel, vel_err = get_predicted_vels (vel_file_list)
    vel0, vel0_n, vel_corr = get_RAVE_vels_corr (corr_file)
    length = len(vel)
    vel0, vel0_n, vel_corr = vel0[0:length], vel0_n[0:length], vel_corr[0:length]
    vel1 = vel + vel_corr
    vel1_err = vel_err
    indeces = np.argwhere(~np.isnan(vel1))
    vel0, vel0_n, vel1 = vel0[indeces], vel0_n[indeces], vel1[indeces]
    print vel0
    print vel1
    plot_stuff (vel0, vel0_n, vel1)
