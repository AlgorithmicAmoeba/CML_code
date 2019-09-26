import pandas
import numpy
import matplotlib.pyplot as plt
from model import Model


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
    Fg_out *= 0

    Cn_in = 0.625 / 60  # (g/L) / (g/mol) = mol/L
    Fn_in = 0.625 / 1000 / Cn_in / 60  # (mg/h) / (mg/g) / (mol/L) / (g/mol) = L/h

    Fb_in = 0.0005  # L/h
    Cb_in = 10  # mol/L

    Fm_in = 0
    Fout = Fg_in + Fn_in + Fb_in + Fm_in

    Tamb = 25
    Q = 5 / 9

    return Fg_in, Cg_in, Fco_in, Cco_in, Fo_in, Co_in, Fg_out, Cn_in, Fn_in, Fb_in, Cb_in, Fm_in, Fout, Tamb, Q


def CgFg(t):
    return numpy.interp(t, glucose['Time'], glucose['Glucose dosing (g/h)'])  # g/h


conc = pandas.read_csv("data/run_7_conc.csv")
glucose = pandas.read_csv("data/run_7_glucose.csv")

print(inputs(0))

ts = numpy.linspace(0, list(glucose['Time'])[-1], 1000)
dt = ts[1]
# Biomass C H_1.8 O_0.5 N_0.2 => 24.6 g/mol
#     Ng, Nx, Nfa, Ne, Nco, No, Nn, Nb, V, Vg
X0 = [0,  4.6/24.6,  0,   0,  0,   0,  0,  0, 4.1,  0.8*2.1, 1.077, 0.1]

Xs = [X0]
m = Model()

for t in ts[1:]:
    dX = numpy.array(m.DEs(Xs[-1], t, inputs))

    X = numpy.array(Xs[-1]) + dX*dt
    Xs.append(list(X))

Xs = numpy.array(Xs)

Vs = Xs[:, 10]
Cgs = Xs[:, 0] * 180 / Vs
Cfas = Xs[:, 2] * 116 / Vs
Ces = Xs[:, 3] * 46 / Vs
Czs = Xs[:, 8] / Vs
Cys = Xs[:, 9] / Vs

plt.figure(figsize=(20, 20))
plt.subplot(2, 2, 1)
plt.plot(ts, Cgs)
plt.plot(conc['Time'], conc['Glucose'], '.')
plt.title("Glucose")

plt.subplot(2, 2, 2)
plt.plot(ts, Cfas)
plt.plot(conc['Time'], conc['Fumaric'], '.')
plt.title("Fumaric")

plt.subplot(2, 2, 3)
plt.plot(ts, Ces)
plt.plot(conc['Time'], conc['Ethanol'], '.')
plt.title("Ethanol")

plt.subplot(2, 2, 4)
plt.plot(ts, Czs, label="Z")
plt.plot(ts, Cys, label="Y")
plt.title("Enzyme")
plt.legend()

plt.show()



