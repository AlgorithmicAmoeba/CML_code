import pandas
import numpy
import matplotlib.pyplot as plt
from model import Model
import StateEstimator
import inputters
import tqdm


inputs = inputters.FakeInputs("data/run_9_glucose.csv")
concentration = pandas.read_csv("data/run_9_conc.csv")

ts = numpy.linspace(0, 200, 200)
# ts = numpy.linspace(0, 30, 30)
# Biomass C H_1.8 O_0.5 N_0.2 => 24.6 g/mol
#     Ng, Nx, Nfa, Ne, Nco, No, Nn, Na, Nb, Nz, Ny, V, Vg
X0 = [0, 4.6/24.6, 0, 0, 0, 0, 0, 1e-5, 0, 5.1, 1.2, 1.077, 0.1]

m = Model(X0, inputs)
Xs = [X0]

# State estimation
t_predict = 1
se = StateEstimator.StateEstimator(X0, inputs, t_predict)

ts_meas = concentration['Time']
Cg_meas, Cfa_meas, Ce_meas = concentration['Glucose'], concentration['Fumaric'], concentration['Ethanol']

t_old_meas = 0
ind_next_measure = 1
t_next_meas = ts_meas[ind_next_measure]

for ti in tqdm.tqdm(ts[1:]):
    # m.step(ti)
    se.step(ts[1])
    if ti > t_next_meas:
        z = [Ci[ind_next_measure] for Ci in [Cg_meas/180, Cfa_meas/116, Ce_meas/46]]
        se.update(z)

        t_old_meas = t_next_meas
        ind_next_measure += 1
        t_next_meas = ts_meas[ind_next_measure]

Xs = se.get_Xs()

Vs = Xs[:, 11]
Cgs = Xs[:, 0] * 180 / Vs
Cfas = Xs[:, 2] * 116 / Vs
Ces = Xs[:, 3] * 46 / Vs
Czs = Xs[:, 9] / Vs
Cys = Xs[:, 10] / Vs
# pH = Xs[:, 13]

Pgs = se.get_Ps()[:, 0] * 180 / Vs

plt.figure(figsize=(20, 20))
plt.subplot(2, 2, 1)
plt.plot(ts, Cgs + Pgs)
plt.plot(ts, Cgs - Pgs)
plt.plot(ts_meas, Cg_meas, '.')
plt.title("Glucose")

plt.subplot(2, 2, 2)
plt.plot(ts, Cfas)
plt.plot(ts_meas, Cfa_meas, '.')
plt.title("Fumaric")

plt.subplot(2, 2, 3)
plt.plot(ts, Ces)
plt.plot(ts_meas, Ce_meas, '.')
plt.title("Ethanol")

plt.subplot(2, 2, 4)
plt.plot(ts, Czs, label="Z")
plt.plot(ts, Cys, label="Y")
plt.title("Enzyme")
plt.legend()

plt.show()

# plt.plot(ts, pH)
# plt.show()


