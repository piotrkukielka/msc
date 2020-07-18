import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
import re
from sklearn.metrics import mean_squared_error
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
    return pd.read_csv(path, sep=" ")


if __name__ == "__main__":
    df_simulations = load_measurepoints("results/")
    df = load_real_data("../real.csv")
    # plt.plot(df_simulations)
    # plt.show()
    # plt.plot(df_real["distance"], df_real["DIC"], '.')
    # plt.show()
    errors = []
    times = []
    distances = []

    for column in df_simulations.iteritems():
        col_name = column[0][0]
        time, distance = re.findall('\d+\.+\d*', col_name)
        times.append(time)
        distances.append(distance)
        errors.append(mean_squared_error(df["DIC"], column[1].values))

    results = pd.DataFrame({'x': distances,
                            't': times,
                            'mse': errors})

    fig = plt.figure(figsize=(16, 9))
    ax = Axes3D(fig)
    surf = ax.plot_trisurf(results['x'], results['t'], results['mse'], cmap=cm.coolwarm)
    # TODO: check x and y
    cbar = fig.colorbar(surf, shrink=0.5, aspect=10)
    cbar.ax.get_yaxis().labelpad = 20
    cbar.ax.set_ylabel('MSE', rotation=270)
    plt.title('Sensitivity analysis')
    plt.xlabel('x [km]')
    plt.ylabel('t [h]')
    plt.show()
