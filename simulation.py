import pandas
import numpy
import matplotlib.pyplot as plt
from model import Model
import StateEstimator
import inputters
import stateUpdaters
import tqdm


inputs = inputters.FakeInputs("data/run_9_glucose.csv")
su = stateUpdaters.FakeStateUpdate("data/run_9_conc.csv")

# ts = numpy.linspace(0, 200, 200)
ts = numpy.linspace(0, 30, 60)
# Biomass C H_1.8 O_0.5 N_0.2 => 24.6 g/mol
#     Ng, Nx, Nfa, Ne, Nco, No, Nn, Na, Nb, Nz, Ny, V, Vg
X0 = [0, 4.6/24.6, 0, 0, 0, 0, 0, 1e-5, 0, 5.1, 1.2, 1.077, 0.1]

m = Model(X0, inputs)
Xs = [X0]

# State estimation
t_predict = 1
se = StateEstimator.StateEstimator(X0, inputs, t_predict)

for ti in tqdm.tqdm(ts[1:]):
    # m.step(ti)
    se.step(ts[1])
    su.step(ts[1])
    if su.update_ready():
        z = su.get_update()
        se.update(z)

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
plt.plot(su.ts_meas, su.Cg_meas, '.')
plt.title("Glucose")

plt.subplot(2, 2, 2)
plt.plot(ts, Cfas)
plt.plot(su.ts_meas, su.Cfa_meas, '.')
plt.title("Fumaric")

plt.subplot(2, 2, 3)
plt.plot(ts, Ces)
plt.plot(su.ts_meas, su.Ce_meas, '.')
plt.title("Ethanol")

plt.subplot(2, 2, 4)
plt.plot(ts, Czs, label="Z")
plt.plot(ts, Cys, label="Y")
plt.title("Enzyme")
plt.legend()

plt.show()

# plt.plot(ts, pH)
# plt.show()


