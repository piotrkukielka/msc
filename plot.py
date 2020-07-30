import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
import re
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D


def load_measurepoints(directory):
    dfs = {}
    os.chdir(directory)
    for filename in glob.glob("*measurepoints"):
        dfs[filename] = pd.read_csv(filename, header=None)
    df = pd.concat(dfs, axis=1)
    return df


def load_real_data(path):
    return pd.read_csv(path, sep=",")


def plot_mse():
    df_simulations = load_measurepoints("results/")
    df = load_real_data("../real_short.csv")
    errors = []
    adv_coeffs = []
    diff_coeffs = []

    for column in df_simulations.iteritems():
        col_name = column[0][0]
        adv_coeff, diff_coeff = re.findall('\d+\.+\d*', col_name) # TODO: check if it is correct
        adv_coeffs.append(adv_coeff)
        diff_coeffs.append(diff_coeff)
        mse = np.square(np.subtract(df["O2"], column[1].values)).mean()
        errors.append(mse)

    results = pd.DataFrame({'adv_coeffs': adv_coeffs,
                            'diff_coeffs': diff_coeffs,
                            'mse': errors})

    fig = plt.figure(figsize=(16, 9))
    ax = Axes3D(fig)
    surf = ax.plot_trisurf(results['adv_coeffs'], results['diff_coeffs'], results['mse'], cmap=cm.coolwarm)
    # TODO: check x and y
    cbar = fig.colorbar(surf, shrink=0.5, aspect=10)
    cbar.ax.get_yaxis().labelpad = 20
    cbar.ax.set_ylabel('MSE', rotation=270)
    plt.title('Sensitivity analysis')
    plt.xlabel('adv_coeffs')
    plt.ylabel('diff_coeffs')
    plt.show()


def plot_one():
    df_real = load_real_data("real.csv")
    print(df_real)
    plt.plot(df_real["distance"], df_real["O2"], '.')
    plt.show()


def plot_short():
    df_real = load_real_data("real_short.csv")
    # df_sim = pd.read_csv("results/-1.000000_0.030000measurepoints", names=["sim_O2"])
    df_sim = pd.read_csv("results/-1.100000_0.100000measurepoints", names=["sim_O2"])
    df = pd.concat([df_real, df_sim], axis=1)
    print(df)
    plt.plot(df["time"], df["O2"], '.', label='pomiary')
    plt.plot(df["time"], df["sim_O2"], 'x', label='symulacja')
    plt.xlabel('czas od startu pomiarów [h]')
    plt.ylabel('stężenie tlenu [g/m3]')
    plt.legend()
    plt.show()


if __name__ == "__main__":
    # plot_mse()
    # plot_one()
    plot_short()
