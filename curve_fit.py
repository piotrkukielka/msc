import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np


def func(x, a, b, c, d):
    res = a * np.sin(b*x + c) + d
    return res


if __name__ == '__main__':
    plt.rcParams.update({'font.size': 14})
    # x_data = np.array([
    #     0,
    #     114.999999997672,
    #     220.000000004657,
    #     339.999999997672,
    #     439.999999998836,
    #     495,
    #     540,
    #     810,
    #     1775.00000000233,
    #     1850.00000000582,
    #     1910.00000000233,
    #     2155.00000000466,
    #     2229.99999999767,
    #     2295,
    #     2509.99999999884,
    #     2785.00000000466,
    #     3115.00000000116,
    #     3205.00000000116,
    #     3354.99999999767,
    #     3630.00000000349,
    #     3739.99999999534,
    #     3995.00000000931,
    #     4469.99999999651,
    #     4614.99999999767,
    #     4745.00000000233,
    # ])/60
    #
    # y_data_temp = np.array([
    #     20.9,
    #     21.6,
    #     21.1,
    #     22.2,
    #     22.8,
    #     23.4,
    #     22.4,
    #     24.1,
    #     21.7,
    #     22.3,
    #     22.6,
    #     23.8,
    #     23.6,
    #     23.4,
    #     22.7,
    #     20.7,
    #     21.1,
    #     21.8,
    #     22.7,
    #     23.8,
    #     23.2,
    #     21.5,
    #     20.2,
    #     21.8,
    #     22.6,
    # ])
    x_data = np.array([2509.99999999884,
              2785.00000000466,
              3115.00000000116,
              3354.99999999767,
              3739.99999999534,
              3995.00000000931,
              4469.99999999651])/60  # jezeli przyjmiemy czas 0 pierwszy pomiar w bondarach
    y_data_o2 = np.array([9.17,
                 6.97,
                 8.47,
                 9.9,
                 10.78,
                 8.72,
                 7.83])
    popt, pcov = curve_fit(func,
                           x_data,
                           # y_data_temp,
                           y_data_o2,
                           bounds=([0, 0, -np.inf, -np.inf],
                                   [5, 0.5, np.inf, np.inf])
                           # bounds=([-np.inf, 0, -np.inf, -np.inf],
                           #         [np.inf, 0.01, np.inf, np.inf]) # co2
                           #bounds=([0,0,0,0], [2., 0.01, 5., 10.]) #o2
                           )
    print(popt)
    # plt.plot(x_data, y_data_temp, 'b.', label='dane')
    plt.plot(x_data, y_data_o2, 'b.', label='dane')
    dense_x = np.arange(min(x_data), max(x_data)+1, 1)
    plt.plot(dense_x, func(dense_x, *popt), 'r-', label='dopasowanie')
    # plt.plot(dense_x, func(dense_x, 1.51, 4.35E-3, -1.63, 22.16)) # to dla temp
    # plt.legend(loc='best', ncol=2)
    plt.legend(loc='best', ncol=1)
    # plt.ylabel('Temperatura [°C]')
    plt.ylabel('DO [$\mathrm{g/m^3}$]')
    plt.xlabel('Czas od początku pomiarów [h]')
    plt.show()
    # plt.savefig("right_boundary.pdf",
    #                bbox_inches='tight',
    #                transparent=True,
    #                pad_inches=0)








    # y_data_dic = np.array([0.003130782776296,
    #               0.003219508173005,
    #               0.003164431899948,
    #               0.003175236814912,
    #               0.003098156756933,
    #               0.003157617430061,
    #               0.003150844506482])
