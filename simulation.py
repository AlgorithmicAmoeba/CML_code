import pandas
import numpy
from model import Model
import StateEstimator
import inputters
import stateUpdaters
import tqdm
import plotting


inputs = inputters.FakeInputs("data/run_9_glucose.csv")
su = stateUpdaters.FakeStateUpdate("data/run_9_conc.csv")

ts = numpy.linspace(0, 200, 200)

# Biomass C H_1.8 O_0.5 N_0.2 => 24.6 g/mol
#     Ng, Nx, Nfa, Ne, Nco, No, Nn, Na, Nb, Nz, Ny, V, Vg
X0 = [0, 4.6/24.6, 0, 0, 0, 0, 0, 1e-5, 0, 5.1, 1.2, 1.077, 0.1, 25]

m = Model(X0, inputs, pH_calculations=True)
Xs = [X0]

# State estimation
t_predict = 1
se = StateEstimator.StateEstimator(X0, inputs, t_predict)

for ti in tqdm.tqdm(ts[1:]):
    m.step(ts[1])
    se.step(ts[1])
    su.step(ts[1])
    if su.update_ready():
        z = su.get_update()
        se.update(z)

Xs = se.get_Xs()
xls = pandas.ExcelWriter('results/result.xls')

model_names = ['Ng', 'Nx', 'Nfa', 'Ne', 'Nco', 'No', 'Nn', 'Na', 'Nb', 'Nz', 'Ny', 'V', 'Vg', 'T', 'pH']
model_data = pandas.DataFrame(m.get_data(), index=ts, columns=model_names)
model_data.index.name = 'ts'
model_data.to_excel(xls, 'model')

se_names = [name + add for add in ['', '_cov'] for name in model_names[:-1]]
se_data = pandas.DataFrame(se.get_data(), index=ts, columns=se_names)
se_data.index.name = 'ts'
se_data.to_excel(xls, 'se')

su_names = ['Cg', 'Cfa', 'Ce']
su_data = pandas.DataFrame(su.get_data().values, index=su.get_times(), columns=su_names)
su_data.index.name = 'ts'
su_data.to_excel(xls, 'su')

xls.save()

plotting.plot_all('results/result.xls')
