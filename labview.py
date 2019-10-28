# Contains all functions that labview will call
import inputters


class Labview:
    def __init__(self):
        self.inputs = inputters.LabviewInputs()


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