import pandas
import numpy


class FakeStateUpdate:
    def __init__(self, update_data_file, t=0, backdate=0):
        self.concentration = pandas.read_csv(update_data_file)
        self.t = t
        self.backdate = backdate

        self.ts_meas = self.concentration['Time']
        self.Cg_meas = self.concentration['Glucose']
        self.Cfa_meas = self.concentration['Fumaric']
        self.Ce_meas = self.concentration['Ethanol']

        self.Cis = [self.Cg_meas / 180, self.Cfa_meas / 116, self.Ce_meas / 46]

        self.t_old_meas = t
        self.ind_next_measure = 1
        self.t_next_meas = self.ts_meas[self.ind_next_measure] + self.backdate

    def step(self, dt):
        self.t += dt

    def update_ready(self):
        return self.t > self.t_next_meas

    def get_update(self):
        if not self.update_ready():
            raise ValueError("Can only get update when one is avaliable")

        z = [Ci[self.ind_next_measure] for Ci in self.Cis]

        self.t_old_meas = self.t_next_meas
        if self.ind_next_measure + 1 < len(self.ts_meas):
            self.ind_next_measure += 1
            self.t_next_meas = self.ts_meas[self.ind_next_measure] + self.backdate
        return z

    def get_times(self):
        return self.ts_meas

    def get_data(self):
        return self.concentration[['Glucose', 'Fumaric', 'Ethanol']].values


class LabviewStateUpdate:
    def __init__(self, t=0):
        self.update = False
        self.update_values = [[0]*3]
        self.update_value = None
        self.update_time = None
        self.ts = [0]

    def step(self, dt):
        pass

    def update_ready(self):
        return self.update

    def get_update(self):
        if not self.update_ready():
            raise ValueError("Can only get update when one is avaliable")

        self.ts.append(self.update_time)
        self.update_values.append(self.update_value)
        return self.update_time, self.update_value

    def get_times(self):
        return numpy.array(self.ts)

    def get_data(self):
        return numpy.array(self.update_values)
