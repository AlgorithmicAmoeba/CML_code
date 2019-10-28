# Contains all functions that labview will call
import inputters
import stateUpdaters
import Model
import StateEstimator


class Labview:
    def __init__(self):
        self.inputs = inputters.LabviewInputs()
        self.su = stateUpdaters.LabviewStateUpdate()
        self.X0 = [0, 4.6/24.6, 0, 0, 0, 0, 0, 1e-5, 0, 5.1, 1.2, 1.077, 0.1, 25]

        self.m = Model.Model(self.X0, self.inputs, pH_calculations=True)
        self.Xs = [self.X0]

        # State estimation
        self.t_predict = 1
        self.se = StateEstimator.StateEstimator(self.X0, self.inputs, self.t_predict)


lv = Labview()


def init():
    """
    Initialises the labview interface.
    Called before the while loop in labview
    """
    pass


def finalise():
    """
    Ends the labview interface.
    Called after the while loop in labview
    """
    pass


def update_inputs(t, inputs):
    """
    Passes inputs into the system to the model.
    Parameters
    ----------
    t : float
        Current time

    inputs : array_like
        The values of the inputs
    """
    lv.inputs.update(t, inputs)


def get_model_outputs():
    """
    Passes outputs to labview from the model
    """
    pass