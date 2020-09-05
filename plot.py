import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
import re
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from sklearn.metrics import mean_squared_error


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
        adv_coeff, diff_coeff = re.findall('\d+\.+\d*', col_name)  # TODO: check if it is correct
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
    modeled = "O2"
    # results_file = "results/-1.000000_0.030000measurepoints"
    results_file = 'results/all.do'
    fig = plt.figure()
    ax = plt.subplot(111)
    df_real = load_real_data("real.csv")
    df_sim = pd.read_csv(results_file, names=["sim_"+modeled])
    df = df_real.join(df_sim)
    if modeled=='DIC':
        df = df.drop(0)
    print(df)
    grouped = df.groupby('place')
    colors=cm.rainbow(np.linspace(0, 1, df['place'].nunique()))
    for (place, data), color in zip(grouped, colors):
        plt.plot(data["time"], data[modeled], 'x', color=color, label=place)
        plt.plot(data["time"], data[modeled], '-', color=color)
        plt.plot(data["time"], data["sim_" + modeled], '.', color=color)
    plt.xlabel('time [h]')
    if modeled == 'O2':
        plt.ylabel('DO [g/m3]')
    if modeled == 'DIC':
        plt.ylabel('DIC [mol/l]')
    # ax.legend()
    # ax.legend(loc='upper left', bbox_to_anchor=(0.3, 1.2),
    #           ncol=3, fancybox=True, shadow=True)
    print(mean_squared_error(df[modeled], df['sim_'+modeled]))
    plt.show()


def plot_short():
    df_real = load_real_data("real_short.csv")
    # df_sim = pd.read_csv("results/-1.000000_0.030000measurepoints", names=["sim_DIC"])
    # df_sim = pd.read_csv("results/-0.000000_0.030000measurepoints", names=["sim_O2"])
    df = pd.concat([df_real, df_sim], axis=1)
    df = df[0:7]
    print(df)
    # plt.plot(df["time"], df["DIC"], '.', label='pomiary')
    # plt.plot(df["time"], df["sim_DIC"], 'x', label='symulacja')
    plt.plot(df["time"], df["O2"], '.', label='pomiary')
    plt.plot(df["time"], df["sim_O2"], 'x', label='symulacja')
    plt.xlabel('czas od startu pomiarów [h]')
    # plt.ylabel('stężenie DIC [mol/L]')
    plt.ylabel('stężenie DO [mol/L]')
    plt.legend()
    plt.show()


def slice():
    names = ['3', '6.init', '9.init']
    for name in names:
        df = pd.read_csv('results/'+name,
                         sep=" ", names=['time', 'distance', 'value'])
        df = df[df['distance']==10]
        plt.plot(df['time'], df['value'])
    plt.show()

if __name__ == "__main__":
    # plot_mse()
    # plot_one()
    # plot_short()
    slice()
