import pandas


class FakeStateUpdate:
    def __init__(self, update_data_file, t=0):
        concentration = pandas.read_csv(update_data_file)
        self.t = t

        self.ts_meas = concentration['Time']
        self.Cg_meas = concentration['Glucose']
        self.Cfa_meas = concentration['Fumaric']
        self.Ce_meas = concentration['Ethanol']

        self.Cis = [self.Cg_meas / 180, self.Cfa_meas / 116, self.Ce_meas / 46]

        self.t_old_meas = t
        self.ind_next_measure = 1
        self.t_next_meas = self.ts_meas[self.ind_next_measure]

    def step(self, dt):
        self.t += dt

    def update_ready(self):
        return self.t > self.t_next_meas

    def get_update(self):
        if not self.update_ready():
            raise ValueError("Can only get update when one is avaliable")

        z = [Ci[self.ind_next_measure] for Ci in self.Cis]

        self.t_old_meas = self.t_next_meas
        self.ind_next_measure += 1
        self.t_next_meas = self.ts_meas[self.ind_next_measure]
        return z
