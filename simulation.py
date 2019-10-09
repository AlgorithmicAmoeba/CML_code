import pandas
import numpy
import matplotlib.pyplot as plt
from model import Model
import tqdm
import filterpy.kalman


def inputs(t):

    Cg_in = 314.19206 / 180  # (g/L) / (g/mol) = mol/L
    Fg_in = CgFg(t) / 180 / Cg_in  # (g/h) / (g/mol) / (mol/L) = L/h

    Qco_in = 8.67 / 1000 * 60  # (ml / min) / (ml/L) * (min/h) = L/h
    Fco_in = 87 * Qco_in / 8.314 / 298  # (kPa) * (L/h) / (L*kPa/mol/K) / (K) = mol/h
    Cco_in = 8.7  # mol CO2 / mol total

    Qo_in = 99.92 / 1000 * 60  # (ml / min) / (ml/L) * (min/h) = L/h
    Fo_in = 87 * Qo_in / 8.314 / 298  # (kPa) * (L/h) / (L*kPa/mol/K) / (K) = mol/h
    Co_in = 21  # mol CO2 / mol total

    Fg_out = Fco_in + Fo_in

    Cn_in = 0.625*10 / 60  # (g/L) / (g/mol) = mol/L
    Fn_in = 0.625 / 1000 / Cn_in / 60  # (mg/h) / (mg/g) / (mol/L) / (g/mol) = L/h

    Fb_in = 0.00006  # L/h
    Cb_in = 10  # mol/L

    Fm_in = 0
    F_out = Fg_in + Fn_in + Fb_in + Fm_in

    T_amb = 25
    Q = 5 / 9

    return Fg_in, Cg_in, Fco_in, Cco_in, Fo_in, Co_in, Fg_out, Cn_in, Fn_in, Fb_in, Cb_in, Fm_in, F_out, T_amb, Q


def CgFg(t):
    return numpy.interp(t, glucose['Time'], glucose['Glucose dosing (g/h)'])  # g/h


concentration = pandas.read_csv("data/run_9_conc.csv")
glucose = pandas.read_csv("data/run_9_glucose.csv")

print(inputs(0))

ts = numpy.linspace(0, list(glucose['Time'])[-1], 1000)
dt = ts[1]
# Biomass C H_1.8 O_0.5 N_0.2 => 24.6 g/mol
#     Ng, Nx, Nfa, Ne, Nco, No, Nn, Na, Nb, Nz, Ny, V, Vg
X0 = numpy.array([0, 4.6/24.6, 0, 0, 0, 0, 0, 1e-5, 0, 5.1, 1.2, 1.077, 0.1])

m = Model(X0)
Xs = [m.outputs()]


# State estimation
def hx(x):
    Ng, _, Nfa, Ne, _, _, _, _, _, _, _, V, _ = x
    return Ng/V, Nfa/V, Ne/V


class fx_obj:
    def __init__(self, t=0):
        self.t = t

    def __call__(self, x, dt):
        ts = numpy.linspace(0, dt, 100)
        dt_small = ts[1]
        m_f = Model(x)
        m_f.t = self.t
        for _ in ts:
            m_f.step(inputs, dt_small)
        return m_f.X


fx = fx_obj()
nx = len(X0)
Q = numpy.diag(numpy.array([1.e-06, 1.e-05, 1.e-05, 1.e-06, 1.e-05, 1.e-05, 1.e-05, 1.e-05,
       1.e-05, 1.e-05, 1.e-05, 1.e-05, 1.e-05]))
R = numpy.diag(numpy.full(3, 1e-5))
sigmas = filterpy.kalman.MerweScaledSigmaPoints(nx, 1e-3, 2, 0)
ukf = filterpy.kalman.UnscentedKalmanFilter(nx, 3, 6 * 3600, hx, fx, sigmas)
ukf.x = X0.copy()
ukf.Q = Q
ukf.R = R

ts_meas = concentration['Time']
Cg_meas, Cfa_meas, Ce_meas = concentration['Glucose'], concentration['Fumaric'], concentration['Ethanol']

t_old_meas = 0
ind_next_measure = 1
t_next_meas = ts_meas[ind_next_measure]
for ti in tqdm.tqdm(ts[1:]):
    Xs.append(list(m.step(inputs, dt)))

    if ti > t_next_meas:
        ukf.predict(t_next_meas - t_old_meas)
        print('p', ukf.x)
        z = [Ci[ind_next_measure] for Ci in [Cg_meas/180, Cfa_meas/116, Ce_meas/46]]
        ukf.update(z)
        print('u', ukf.x)
        m.X = ukf.x.copy()

        t_old_meas = t_next_meas
        ind_next_measure += 1
        t_next_meas = ts_meas[ind_next_measure]
        fx.t = ti

Xs = numpy.array(Xs)

Vs = Xs[:, 11]
Cgs = Xs[:, 0] * 180 / Vs
Cfas = Xs[:, 2] * 116 / Vs
Ces = Xs[:, 3] * 46 / Vs
Czs = Xs[:, 9] / Vs
Cys = Xs[:, 10] / Vs
pH = Xs[:, 13]

plt.figure(figsize=(20, 20))
plt.subplot(2, 2, 1)
plt.plot(ts, Cgs)
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

plt.plot(ts, pH)
plt.show()


