# Contains all functions that labview will call

def init():
    """
    Initialises the labview interface.
    Called before the while loop in labview
    """
    pass

def finilise():
    """
    Ends the labview interface.
    Called after the while loop in labview
    """
    pass


def set_model_inputs(inputs):
    """
    Passes inputs into the system to the model.
    Parameters
    ----------
    inputs : array_like
        The values of the inputs
    """
    pass


def get_model_outputs():
    """
    Passes outputs to labview from the model
    """
    pass