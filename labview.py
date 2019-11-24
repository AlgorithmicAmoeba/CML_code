# Contains all functions that labview will call
import inputters
import stateUpdaters
import Model
import StateEstimator
import scipy.stats
import matplotlib.pyplot as plt
import plotting


class Labview:
    """Contains the required objects to store state for labview functions

    Attributes
    ----------
    t : float
        The current time

    ts : array_like
        List of times

    inputs : inputters.LabviewInputs
        Input object that stores and retrieves inputs

    su : stateUpdaters.LabviewStateUpdate
        State update object that stores and retrieves HPLC input data

    X0 : array_like
        The initial state

    m : Model.Model
        The system model

    Xs : array_like
        List of previous states

    t_predict : float
        The period between state predictions

    se : StateEstimator.StateEstimator
        State estimator object

    live_plot : bool
        If `True` then a live plot of the run is shown

    """
    def __init__(self):
        self.t = 0
        self.ts = [self.t]

        self.inputs = inputters.LabviewInputs()
        self.su = stateUpdaters.LabviewStateUpdate()
        self.X0 = [0, 4.6/24.6, 0, 0, 0, 0, 0, 1e-5, 0, 5.1, 1.2, 1.077, 0.1, 35]

        self.m = Model.Model(self.X0, self.inputs, pH_calculations=True)
        self.Xs = [self.X0]

        # State estimation
        self.t_predict = 0.9/3600
        self.se = StateEstimator.StateEstimator(self.X0, self.inputs, self.t_predict)

        # Plotting
        self.live_plot = True


lv = Labview()

if lv.live_plot:
    plt.figure(figsize=(20, 20))
    plt.ion()


def init():
    """ Initialises the labview interface.
    Called before the while loop in labview
    """
    pass


def finalise():
    """Ends the labview interface.
    Called after the while loop in labview
    """
    pass


def update_inputs(t, inputs):
    """ Passes inputs into the system to the model.

    Parameters
    ----------
    t : float
        Current time

    inputs : array_like
        The values of the inputs
    """
    lv.inputs.update(t, inputs)


def step(t):
    """Steps the labview object through time

    Parameters
    ----------
    t : float
        Current time
    """
    dt = t - lv.t
    lv.t = t
    lv.ts.append(t)
    lv.m.step(dt)
    lv.se.step(dt)
    lv.su.step(dt)
    if lv.su.update_ready():
        t_u, z = lv.su.get_update()
        lv.se.update(z, t-t_u)
        lv.su.update = False
        
    if lv.live_plot:
        plotting.plot_live(lv.ts, lv.m, lv.se, lv.su)


def update_state(t, z):
    """Sends update to the state updater

    Parameters
    ----------
    t : float
        Timestamp of the update

    z : array_like
        Update values
    """
    lv.su.update = True
    lv.su.update_value = [zi/n for zi, n in zip(z, [180, 116, 46])]
    lv.su.update_time = t


def get_glucose_graph(confidence=0.95):
    """ Passes outputs to labview from the model

    Parameters
    ----------
    confidence : float
        The confidence probability for the plots
    """
    model = lv.m.get_data()
    se = lv.se.get_data()
    su = lv.su.get_data().values

    ts = lv.ts
    # Model
    Vs_m = model[:, 11]
    Cgs_m = model[:, 0] * 180 / Vs_m

    # State estimator
    Vs = se[:, 11]
    Cgs = se[:, 0] * 180 / Vs_m

    # Standard deviation multiplier to get the correct confidence interval
    K = scipy.stats.norm.ppf(confidence)
    # SE covs
    Pgs = se[:, 14 + 0] * 180 / Vs * K

    # Measured update values
    ts_meas = lv.su.get_times()
    Cg_meas = su[:, 0][:len(ts_meas)]

    plots = [(ts, Cgs_m), (ts, Cgs + Pgs), (ts, Cgs - Pgs), (ts_meas, Cg_meas)]

    return plots
