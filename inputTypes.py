import numpy
import pandas


class FakeInputs:
    def __init__(self, glucode_data_file):
        self.glucose = pandas.read_csv(glucode_data_file)

    def __call__(self, t):
        Cg_in = 314.19206 / 180  # (g/L) / (g/mol) = mol/L
        Fg_in = self.CgFg(t) / 180 / Cg_in  # (g/h) / (g/mol) / (mol/L) = L/h

        Qco_in = 8.67 / 1000 * 60  # (ml / min) / (ml/L) * (min/h) = L/h
        Fco_in = 87 * Qco_in / 8.314 / 298  # (kPa) * (L/h) / (L*kPa/mol/K) / (K) = mol/h
        Cco_in = 8.7  # mol CO2 / mol total

        Qo_in = 99.92 / 1000 * 60  # (ml / min) / (ml/L) * (min/h) = L/h
        Fo_in = 87 * Qo_in / 8.314 / 298  # (kPa) * (L/h) / (L*kPa/mol/K) / (K) = mol/h
        Co_in = 21  # mol CO2 / mol total

        Fg_out = Fco_in + Fo_in

        Cn_in = 0.625 * 10 / 60  # (g/L) / (g/mol) = mol/L
        Fn_in = 0.625 / 1000 / Cn_in / 60  # (mg/h) / (mg/g) / (mol/L) / (g/mol) = L/h

        Fb_in = 0.00006  # L/h
        Cb_in = 10  # mol/L

        Fm_in = 0
        F_out = Fg_in + Fn_in + Fb_in + Fm_in

        T_amb = 25
        Q = 5 / 9

        return Fg_in, Cg_in, Fco_in, Cco_in, Fo_in, Co_in, Fg_out, Cn_in, Fn_in, Fb_in, Cb_in, Fm_in, F_out, T_amb, Q

    def CgFg(self, t):
        return numpy.interp(t, self.glucose['Time'], self.glucose['Glucose dosing (g/h)'])  # g/h

