import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np


def func(x, a, b, c, d):
    res = a * np.sin(b*x + c) + d
    return res


if __name__ == '__main__':
    x_data = np.array([2509.99999999884,
              2785.00000000466,
              3115.00000000116,
              3354.99999999767,
              3739.99999999534,
              3995.00000000931,
              4469.99999999651])  # jezeli przyjmiemy czas 0 pierwszy pomiar w bondarach
    y_data_o2 = np.array([9.17,
                 6.97,
                 8.47,
                 9.9,
                 10.78,
                 8.72,
                 7.83])
    y_data_dic = np.array([0.003130782776296,
                  0.003219508173005,
                  0.003164431899948,
                  0.003175236814912,
                  0.003098156756933,
                  0.003157617430061,
                  0.003150844506482])
    popt, pcov = curve_fit(func,
                           x_data,
                           y_data_dic,
                           bounds=([-np.inf, 0, -np.inf, -np.inf],
                                   [np.inf, 0.01, np.inf, np.inf]) # co2
                           #bounds=([0,0,0,0], [2., 0.01, 5., 10.]) #o2
                           )
    print(popt)
    plt.plot(x_data, y_data_dic, 'b.', label='data')
    dense_x = np.arange(min(x_data), max(x_data), 1)
    plt.plot(dense_x, func(dense_x, *popt), 'r-', label='fit')
    plt.legend()
    plt.show()
