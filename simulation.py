import pandas
import numpy
from Model import Model
import StateEstimator
import inputters
import stateUpdaters
import tqdm
import plotting
import matplotlib.pyplot as plt

backdate = 0
inputs = inputters.FakeInputsReal() # ("data/run_9_glucose.csv")
su = stateUpdaters.FakeStateUpdate("data/run_9_conc.csv", backdate=backdate)

ts = numpy.linspace(0, 230, 1000)

# Biomass C H_1.8 O_0.5 N_0.2 => 24.6 g/mol
#     Ng, Nx, Nfa, Ne, Nco, No, Nn, Na, Nb, Nez, Nfaz, Nezfa, V, Vg, T
X0 = [3.1/180, 1e-3/24.6, 0, 0, 0, 0, 2/60, 1e-5, 0, 0, 0, 0, 1.077, 0.1, 25]

m = Model(X0, inputs, pH_calculations=True)
Xs = [X0]

# State estimation
t_predict = 1
se = StateEstimator.StateEstimator(X0, inputs, t_predict)

live_plot = False

if live_plot:
    plt.figure(figsize=(20, 20))
    plt.ion()

reset = True
for ti in tqdm.tqdm(ts[1:]):
    m.step(ts[1])
    # se.step(ts[1])
    V = m.X[-3]
    # if ti > 20:
    #     m.X[[0, 1, 2, 3, 4, 5, 6]] = numpy.array([0.05/180, 0.6/24.6, 0, 1, 1, 1, 0.001])*V
    # else:
    #     m.X[[0, 1, 2, 3, 4, 5, 6]] = numpy.array([0.5 / 180, 0.6 / 24.6, 0, 1, 1, 1, 0.001]) * V
    if reset and ti > 26:
        reset = False
        m.X[[0, 2, 3, 4, 5, 6]] = 0
        # m.X[6] = 0

    su.step(ts[1])
    if su.update_ready():
        z = su.get_update()
        se.update(z, ti-backdate if backdate else numpy.nan)

    if live_plot:
        plotting.plot_live(ts, m, se, su)

if live_plot:
    plt.ioff()

Xs = se.get_Xs()
xls = pandas.ExcelWriter('results/result.xlsx', engine='xlsxwriter')

model_names = ['Ng', 'Nx', 'Nfa', 'Ne', 'Nco', 'No', 'Nn', 'Na', 'Nb', 'Nez', 'Nfaz', 'Nezfa', 'V', 'Vg', 'T', 'pH']
model_data = pandas.DataFrame(m.get_data(), index=ts, columns=model_names)
model_data.index.name = 'ts'
model_data.to_excel(xls, 'model')
#
# se_names = [name + add for add in ['', '_cov'] for name in model_names[:-1]]
# se_data = pandas.DataFrame(se.get_data(), index=ts, columns=se_names)
# se_data.index.name = 'ts'
# se_data.to_excel(xls, 'se')
#
su_names = ['Cg', 'Cfa', 'Ce']
su_data = pandas.DataFrame(su.get_data(), index=su.get_times(), columns=su_names)
su_data.index.name = 'ts'
su_data.to_excel(xls, 'su')
#
xls.save()
#
plotting.plot_model('results/result.xlsx')
