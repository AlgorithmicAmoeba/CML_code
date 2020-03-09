import numpy
import pandas


class FakeInputs:
    """Creates fake inputs for the glucose feed from past data

    Parameters
    -----------
    glucose_data_file : string
        The name of the csv file that contains the glucose data

    Attributes
    -----------
    glucose : pandas.Dataframe
        An object containing the info from the glucose file
    """
    def __init__(self, glucose_data_file):
        self.glucose = pandas.read_csv(glucose_data_file)

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
        """Interpolates the value from the glucose file

        Parameters
        ----------
        t : float
            The value of time at which the input should be looked up
        """
        return numpy.interp(t, self.glucose['Time'], self.glucose['Glucose dosing (g/h)'])  # g/h


class FakeInputsReal:
    """Creates fake inputs for the glucose feed from past data
    """

    def __call__(self, t):
        t_batch = 30

        Cn_in = 0.625 * 10 / 60  # (g/L) / (g/mol) = mol/L
        if t < t_batch + 66:
            CgFg = 0.141612826257827
        elif t < t_batch + 66 + 80:
            CgFg = 0.21241923938674
        else:
            CgFg = 0.21241923938674 * 2
        Cg_in = 314.19206 / 180  # (g/L) / (g/mol) = mol/L
        Cb_in = 10  # mol/L
        Fm_in = 0

        if t < t_batch:
            Fg_in = 0
            Fn_in = 0
            Fb_in = 0
            F_out = 0
        else:
            Fg_in = CgFg / 180 / Cg_in  # (g/h) / (g/mol) / (mol/L) = L/h
            Fn_in = 0.625 / 1000 / Cn_in / 60  # (mg/h) / (mg/g) / (mol/L) / (g/mol) = L/h
            Fb_in = 0.00006  # L/h
            F_out = Fg_in + Fn_in + Fb_in + Fm_in

        Qco_in = 8.67 / 1000 * 60  # (ml / min) / (ml/L) * (min/h) = L/h
        Fco_in = 87 * Qco_in / 8.314 / 298  # (kPa) * (L/h) / (L*kPa/mol/K) / (K) = mol/h
        Cco_in = 8.7  # mol CO2 / mol total

        Qo_in = 99.92 / 1000 * 60  # (ml / min) / (ml/L) * (min/h) = L/h
        Fo_in = 87 * Qo_in / 8.314 / 298  # (kPa) * (L/h) / (L*kPa/mol/K) / (K) = mol/h
        Co_in = 21  # mol CO2 / mol total

        Fg_out = Fco_in + Fo_in

        T_amb = 25
        Q = 5 / 9

        return Fg_in, Cg_in, Fco_in, Cco_in, Fo_in, Co_in, Fg_out, Cn_in, Fn_in, Fb_in, Cb_in, Fm_in, F_out, T_amb, Q

    def CgFg(self, t):
        """Interpolates the value from the glucose file

        Parameters
        ----------
        t : float
            The value of time at which the input should be looked up
        """
        return numpy.interp(t, self.glucose['Time'], self.glucose['Glucose dosing (g/h)'])  # g/h


class LabviewInputs:
    """Stores and looks up input values from Labview

    Attributes
    -----------
    ts : array_like
        Stores the time stamp information about the inputs

    inputs : array_like
        Stores the inputs

    Cg_in : float, constant
        The glucose feed concentration

    G_rpm_to_ml_min : float, constant
        The conversion factor between rpm and ml/min for the glucose pump

    Cco_in : float, constant
        The percentage CO2 in the CO2 feed

    Co_in : float, constant
        The percentage O2 in the O2 feed

    Cn_in : float, constant
        The concentration urea in the nitrogen feed

    N_rpm_to_ml_min : float, constant
        The conversion factor between rpm and ml/min for the nitrogen pump

    B_rpm_to_ml_min : float, constant
        The conversion factor between rpm and ml/min for the sodium hydroxide pump

    Cb_in : float, constant
        The concentration of sodium hydroxide in the base feed

    M_rpm_to_ml_min : float, constant
        The conversion factor between rpm and ml/min for the mineral pump

    T_amb : float
        The ambient room temperature

    Q_fact : float, constant
        A multiplier for the heater gain
    """
    def __init__(self):
        self.ts = []
        self.inputs = []
        self.Cg_in = 314.19206 / 180  # (g/L) / (g/mol) = mol/L
        self.G_rpm_to_ml_min = 0.02117909  # (ml/min) / (rpm)

        self.Cco_in = 8.7  # mol CO2 / mol total
        self.Co_in = 21  # mol CO2 / mol total

        self.Cn_in = 0.625 * 10 / 60  # (g/L) / (g/mol) = mol/L
        self.N_rpm_to_ml_min = self.G_rpm_to_ml_min

        self.B_rpm_to_ml_min = self.G_rpm_to_ml_min
        self.Cb_in = 10  # mol/L

        self.M_rpm_to_ml_min = 0.01878002

        self.T_amb = 25
        self.Q_fact = 1

    def update(self, t, data):
        """Update the current inputss

        Parameters
        ----------
        t : float
            Current time
        data : array_like
            Current inputs
        """
        self.ts.append(t)
        self.inputs.append(data)

    def __call__(self, t):
        offsets = [0.004422699015892, 0.004053967439397, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0, 0]
        spans = [0.0004717700792543, 9.925512172352e-5, 4e-5, 0.000542, 4e-5, 0.00032, 0.00064, 0.0008, 1, 1]

        index = min(numpy.searchsorted(self.ts, t), len(self.ts)-1)
        ins = [(in_i - offset)/span for in_i, offset, span in zip(self.inputs[index], offsets, spans)]
        CO2_ml_min, O2_ml_min, _, B_rpm, _, M_rpm, G_rpm, N_rpm, B_on_off, Q_on_off = ins

        Cg_in = self.Cg_in
        Fg_in = G_rpm * self.G_rpm_to_ml_min / 1000 * 60  # (ml/min) / (L/ml) * (min/h) = L/h

        Qco_in = CO2_ml_min / 1000 * 60  # (ml / min) / (ml/L) * (min/h) = L/h
        Fco_in = 87 * Qco_in / 8.314 / 298  # (kPa) * (L/h) / (L*kPa/mol/K) / (K) = mol/h
        Cco_in = self.Cco_in

        Qo_in = O2_ml_min / 1000 * 60  # (ml / min) / (ml/L) * (min/h) = L/h
        Fo_in = 87 * Qo_in / 8.314 / 298  # (kPa) * (L/h) / (L*kPa/mol/K) / (K) = mol/h
        Co_in = self.Co_in

        Fg_out = Fco_in + Fo_in

        Cn_in = self.Cn_in
        Fn_in = N_rpm * self.N_rpm_to_ml_min / 1000 * 60  # (ml/min) / (L/ml) * (min/h) = L/h

        Fb_in = B_on_off * B_rpm * self.B_rpm_to_ml_min / 1000 * 60  # (ml/min) / (L/ml) * (min/h) = L/h
        Cb_in = self.Cb_in

        Fm_in = M_rpm * self.M_rpm_to_ml_min
        F_out = Fg_in + Fn_in + Fb_in + Fm_in

        T_amb = self.T_amb
        Q = Q_on_off * self.Q_fact

        return Fg_in, Cg_in, Fco_in, Cco_in, Fo_in, Co_in, Fg_out, Cn_in, Fn_in, Fb_in, Cb_in, Fm_in, F_out, T_amb, Q

    def get_data(self):
        """Get all the input data

        Returns
        -------
        out : array_like
            All the input data
        """
        ts = numpy.array(self.ts)
        inputs = numpy.array(self.inputs)
        out = numpy.concatenate([ts, inputs])
        return out
