import matplotlib.pyplot as plt
import pandas
import scipy.stats


def plot_all(file_name, confidence=0.95, show=True):
    xls = pandas.ExcelFile(file_name)
    model = pandas.read_excel(xls, 'model')
    se = pandas.read_excel(xls, 'se')
    su = pandas.read_excel(xls, 'su')

    # SE means
    ts = se['ts']
    Vs = se['V']
    Cgs = se['Ng'] * 180 / Vs
    Cfas = se['Nfa'] * 116 / Vs
    Ces = se['Ne'] * 46 / Vs
    Czs = se['Nz'] / Vs
    Cys = se['Ny'] / Vs

    # Standard deviation multiplier to get the correct confidence interval
    K = scipy.stats.norm.ppf(confidence)
    # SE covs
    Pgs = se['Ng_cov'] * 180 / Vs * K
    Pfas = se['Nfa_cov'] * 116 / Vs * K
    Pes = se['Ne_cov'] * 46 / Vs * K
    Pzs = se['Nz_cov'] / Vs * K
    Pys = se['Ny_cov'] / Vs * K

    # Measured update values
    ts_meas = su['ts']
    Cg_meas = su['Cg']
    Cfa_meas = su['Cfa']
    Ce_meas = su['Ce']

    plt.figure(figsize=(20, 20))
    plt.subplot(2, 2, 1)
    plt.plot(ts, Cgs + Pgs)
    plt.plot(ts, Cgs - Pgs)
    plt.plot(ts_meas, Cg_meas, '.')
    plt.title("Glucose")

    plt.subplot(2, 2, 2)
    plt.plot(ts, Cfas + Pfas)
    plt.plot(ts, Cfas - Pfas)
    plt.plot(ts_meas, Cfa_meas, '.')
    plt.title("Fumaric")

    plt.subplot(2, 2, 3)
    plt.plot(ts, Ces + Pes)
    plt.plot(ts, Ces - Pes)
    plt.plot(ts_meas, Ce_meas, '.')
    plt.title("Ethanol")

    plt.subplot(2, 2, 4)
    plt.plot(ts, Czs + Pzs, label="Z+")
    plt.plot(ts, Czs - Pzs, label="Z-")

    plt.plot(ts, Cys + Pys, label="Y+")
    plt.plot(ts, Cys - Pys, label="Y-")
    plt.title("Enzyme")
    plt.legend()

    if show:
        plt.show()

    # plt.plot(ts, pH)
    # plt.show()


# plot_all('data/result.xls')
