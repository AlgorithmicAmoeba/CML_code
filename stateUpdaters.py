import pandas
import numpy


class FakeStateUpdate:
    """Creates fake state updates from past data

    Parameters
    -----------
    update_data_file : string
        The name of the csv file that contains the data

    t : float, optional
        Initial time.
        Defaults to zero

    backdate : float, optional
        Sets a time to back date samples.
        Simulates the fact that HPLC readings are not instant.
        Defaults to zero

    Attributes
    -----------
    concentration : pandas.Dataframe
        An object containing the info from the data file

    t : float
        Current time.

    backdate : float
        Sets a time to back date samples

    ts_meas, Cg_meas, Cfa_meas, Ce_meas : array_like
        List of HPLC readings and time stamps

    Cis : array_like
        Adjusted HPLC readings
    """
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
        """Steps the updater through time
        Parameters
        ----------
        dt : float
            Time since the previous step
        """
        self.t += dt

    def update_ready(self):
        """Returns `True` if there is an update ready for the state estimator
        """
        return self.t > self.t_next_meas

    def get_update(self):
        """Get the state update for the state estimator

        Returns
        --------
        z : array_like
            The update
        """
        if not self.update_ready():
            raise ValueError("Can only get update when one is avaliable")

        z = [Ci[self.ind_next_measure] for Ci in self.Cis]

        self.t_old_meas = self.t_next_meas
        if self.ind_next_measure + 1 < len(self.ts_meas):
            self.ind_next_measure += 1
            self.t_next_meas = self.ts_meas[self.ind_next_measure] + self.backdate
        return z

    def get_times(self):
        """Gets the times of the state updates
        """
        return self.ts_meas

    def get_data(self):
        """Gets all data from the object"""
        return self.concentration[['Glucose', 'Fumaric', 'Ethanol']].values


class LabviewStateUpdate:
    """Handles the state updates from Labview

    Attributes
    ----------
    update: bool
        `True` if there is an update ready for the state estimator

    update_values : array_like
        List of previous update values

    update_value : array_like
        List of current update value

    update_time : float
        Time at which the readings were taken

    ts : array_like
        List of previous times
    """
    def __init__(self):
        self.update = False
        self.update_values = [[0]*3]
        self.update_value = None
        self.update_time = None
        self.ts = [0]

    def step(self, dt):
        """Steps the updater through time"""
        pass

    def update_ready(self):
        """Returns `True` if there is an update ready for the state estimator
        """
        return self.update

    def get_update(self):
        """Get the state update for the state estimator

        Returns
        --------
        update_time : float
            Time at which the readings were taken

        update_value : array_like
            The update
        """
        if not self.update_ready():
            raise ValueError("Can only get update when one is avaliable")

        self.ts.append(self.update_time)
        self.update_values.append([zi*n for zi, n in zip(self.update_value, [180, 116, 46])])
        return self.update_time, self.update_value

    def get_times(self):
        """Gets the times of the state updates
        """
        return numpy.array(self.ts)

    def get_data(self):
        """Gets all data from the object"""
        return numpy.array(self.update_values)
