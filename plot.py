import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
import re
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from sklearn.metrics import mean_squared_error
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D


def load_measurepoints(directory):
    dfs = {}
    os.chdir(directory)
    for filename in glob.glob("*measurepoints"):
        dfs[filename] = pd.read_csv(filename, header=None)
    df = pd.concat(dfs, axis=1)
    return df

def plot_3d():
    df = pd.read_csv('real_results/old_preparameters/peclet/pe2.4', names=['t', 'x', 'z'], sep=' ')
    fig = plt.figure(figsize=(16, 10))
    ax = Axes3D(fig)
    surf = ax.plot_trisurf(df['t'], df['x'], df['z'], antialiased=False, cmap=cm.magma)
    cbar = fig.colorbar(surf, shrink=0.5, aspect=10, pad=-0.2)
    cbar.ax.set_ylabel('stężenie', rotation=270)
    plt.xlabel('czas')
    plt.ylabel('odległość')
    ax.set_zlabel('stężenie')
    ax.view_init(elev=30., azim=250)
    # ax.set_ylim(bottom=0)
    ax.set_xlim(1, 1.2)
    ax.set_ylim(3, 3.15)
    # plt.show()
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    plt.savefig("2_4.pdf",
                bbox_inches='tight',
                transparent=True,
                pad_inches=0)


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


def export_legend(legend, filename="legend.pdf"):
    fig  = legend.figure
    fig.canvas.draw()
    bbox = legend.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(filename, dpi="figure", bbox_inches=bbox)


def plot_one(results_file):
    plt.rcParams.update({'font.size': 11})
    modeled = "DIC"
    # results_file = "results/-1.000000_0.030000measurepoints"
    # results_file = 'results/normal.030000measurepoints'
    if modeled is "O2":
        error_multiplier = 0.03
    if modeled is "DIC":
        error_multiplier = 0.05
    fig = plt.figure()
    ax = plt.subplot(111)
    df_real = load_real_data("real.csv")
    # df_real = df_real.dropna()
    df_sim = pd.read_csv(results_file, names=["sim_" + modeled])
    df = df_real.join(df_sim)
    df = df.dropna()
    print(df)
    grouped = df.groupby('place')
    # colors=cm.Set1(np.linspace(0, 1, df['place'].nunique()))
    colors = cm.tab10.colors[0:df['place'].nunique()]
    patches = []
    for (place, data), color in zip(grouped, colors):
        plt.plot(data["time"], data[modeled], '_', color=color, label=place, markersize=5)
        plt.plot(data["time"], data[modeled], '--', color=color, linewidth=1, alpha=1)
        plt.plot(data["time"], data["sim_" + modeled], '.', color=color, markersize=9)
        plt.errorbar(x=data["time"], y=data[modeled], yerr=data[modeled] * error_multiplier, capsize=5, fmt='none')
        patches.append(mpatches.Patch(color=color, label=place))
    plt.xlabel('czas od początku symulacji [h]')
    if modeled == 'O2':
        plt.ylabel('DO [$\mathrm{g/m^3}$]')
    if modeled == 'DIC':
        plt.ylabel('DIC [mol/l]')
    legend_patches = plt.legend(handles=patches,
                                loc='best',
                                ncol=3)#, fancybox=True, shadow=True)
    markers = [Line2D([], [], color='black', lw=0, label='Pomiar', marker='_'),
               Line2D([], [], color='black', lw=0, label='Symulacja', marker='x')]
    # legend_location = 'best' if modeled=="O2" else (0.01, 0.65)
    legend_location = 'best' if modeled=="O2" else (0.01, 0.88) # for dic with r5 xD
    legend_marker = plt.legend(handles=markers,
                               loc=legend_location,
                               ncol=1)#, fancybox=True, shadow=True)
    # possible on one legend if needed
    print(mean_squared_error(df[modeled], df['sim_' + modeled]))
    plt.gca().add_artist(legend_patches)
    plt.gca().add_artist(legend_marker)
    fig.set_size_inches(10, 6)
    # plt.show()
    plt.savefig(results_file+".pdf",
                   bbox_inches='tight',
                   transparent=True,
                   pad_inches=0)
    #
    # export_legend(legend_patches)
    # plt.show()
    # plt.savefig("dic_data_original.pdf",
    #                bbox_inches='tight',
    #                transparent=True,
    #                pad_inches=0)


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


def slice(slice_place, thing, on_what):
    plt.rcParams.update({'font.size': 14})
    if on_what == 'distance':
        x_axis = 'time'
    if on_what == 'time':
        x_axis = 'distance'
    # names = ['3', '6'+thing, '9'+thing, '12'+thing]
    names = ['3', '6', '9', '12']
    for name in names:
        name_solo = name
        name = name + thing
        df = pd.read_csv('results/' + name,
                         sep=" ", names=['time', 'distance', 'value'])
        df = df[df[on_what] == slice_place]
        plt.plot(df[x_axis], df['value'], label=name_solo + ' $\mathrm{g/m^3}$')
    plt.legend()
    plt.ylabel('DO [$\mathrm{g/m^3}$]')
    if on_what == 'distance':
        plt.xlabel('czas od początku symulacji [h]')
        plt.xlim(0, 100)
    if on_what == 'time':
        plt.xlabel('odległość [km]')
    # plt.show()
    plt.savefig('results/' + thing + '_slice_' + str(slice_place) + '_on_' + on_what + '.pdf')
    plt.clf()

def slice_right_bound():
    on_what = "time"
    x_axis = 'distance'
    slice_place = 180
    names = ["-1.140000_0.080000_rightbound_neumann", "-1.140000_0.080000_rightbound_sinus"]
    for name in names:
        name_solo = name
        name = name
        label = 'warunek Dirichleta' if name=="-1.140000_0.080000_rightbound_sinus" else 'warunek von Neumanna'
        df = pd.read_csv('results/' + name,
                         sep=" ", names=['time', 'distance', 'value'])
        df = df[df[on_what] == slice_place]
        plt.plot(df[x_axis], df['value'], label=label)
    plt.legend()
    plt.ylabel('DO [$\mathrm{g/m^3}$]')
    if on_what == 'distance':
        plt.xlabel('czas od początku symulacji [h]')
        plt.xlim(0, 100)
    if on_what == 'time':
        plt.xlabel('odległość [km]')
    # plt.show()
    plt.savefig('results/' + '_slice_' + str(slice_place) + '_on_' + on_what + '.pdf')
    plt.clf()

def plot_photosynthesis():
    plt.rcParams.update({'font.size': 11})
    fig = plt.figure()
    ax = plt.subplot(111)
    df = pd.read_csv("results/photosynthesis.dat", names=["time", "O2"], sep=" ")
    print(df)
    # dla godziny zegarowej
    df = df[df["time"] <= 24]
    # hr_of_day = (df["time"]+5)%24
    # plt.plot(hr_of_day, df["O2"])
    # plt.plot(hr_of_day, [0]*len(df["time"]), '--')
    # plt.xlabel("godzina dnia (zegarowa)")
    # plt.ylabel("produkcja tlenu [$\mathrm{g/m^3}$]")
    # dla czasu symulacji
    plt.plot(df["time"], df["O2"])
    plt.plot(df["time"], [0]*len(df["time"]), '--')
    plt.xlabel('czas [h]')
    plt.ylabel("produkcja tlenu [$\mathrm{g/m^3}$]")
    fig.set_size_inches(10, 6)
    # plt.show()
    plt.savefig("results/photosynthesis.pdf",
                bbox_inches='tight',
                transparent=True,
                pad_inches=0)



if __name__ == "__main__":
    # plot_3d()
    plot_one("results/-1.140000_0.080000measurepoints")
    # plot_one("results/r2_r5_dic")
    # for filename in glob.glob("results/*"):
    #     plot_one(filename)
    # plot_mse()
    # plot_short()

    # _things = ['.bound', '.init']
    # _slice_places = [10, 50, 90]
    # _on_whats = ['distance', 'time']

    # _things = ['_longboth.bound', '.init']
    # _slice_places = [10, 50, 90, 120, 150, 180]
    # _on_whats = ['distance', 'time']
    # for _thing in _things:
    #     for _slice_place in _slice_places:
    #         for _on_what in _on_whats:
    #             slice(_slice_place, _thing, _on_what)

    # slice_right_bound()

    # plot_photosynthesis()
