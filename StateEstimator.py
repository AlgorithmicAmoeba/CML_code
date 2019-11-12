import numpy
import scipy
import filterpy.kalman
from AdjMerweScaledSigmaPoints import MerweScaledSigmaPoints
import Model


class StateEstimator:
    def __init__(self, X0, inputs, t_predict):
        self.inputs = inputs

        self._Xs = [X0]
        self._Ps = [numpy.zeros((len(X0), len(X0)))]
        self._deviations = [[0] * len(X0)]
        self.ts = [0]

        #                           Ng, Nx, Nfa, Ne, Nco, No, Nn, Na, Nb, Nz, Ny, V, Vg, T
        self.Q = numpy.diag(numpy.array([1e-6, 1e-3, 1e-5, 1e-4, 1e-5, 1e-5, 1e-5,
                                         1e-5, 1e-5, 1e-2, 1e-2, 1e-5, 1e-5, 1e-1]))
        self.R = numpy.diag(numpy.array([1e-12, 1e-12, 1e-12]))

        self.fx = self.FXObj(self.inputs)
        self.nx = len(self.Q)

        self.sigmas = MerweScaledSigmaPoints(self.nx, 1e-3, 2, 0, sqrt_method=scipy.linalg.sqrtm)
        self.ukf = filterpy.kalman.UnscentedKalmanFilter(self.nx, 3, 0, self.hx, self.fx, self.sigmas)

        self.ukf.x = X0
        self.ukf.Q = self.Q
        self.ukf.R = self.R

        self.t = 0
        self.t_next_predict = 0
        self.t_predict = t_predict

        self.t_next_predicts = [self.t_next_predict]

    @staticmethod
    def hx(x):
        Ng, _, Nfa, Ne, _, _, _, _, _, _, _, V, _, _ = x
        return Ng/V, Nfa/V, Ne/V

    class FXObj:
        def __init__(self, inputs, t=0):
            self.t = t
            self.inputs = inputs

        def __call__(self, x, dt):
            ts = numpy.linspace(0, dt, int(dt*5 + 2))
            dt_small = ts[1]
            m_f = Model.Model(x, self.inputs)
            m_f.t = self.t
            for _ in ts:
                m_f.step(dt_small)
            return m_f.X

    def step(self, dt):
        self.t += dt

        if self.t > self.t_next_predict:
            self.fx.t = self.t
            self.ukf.predict(self.t_predict)
            self.t_next_predict = self.t + self.t_predict

        self._Xs.append(self.ukf.x)
        self._Ps.append(self.ukf.P)
        self._deviations.append(numpy.sqrt(numpy.diag(self.ukf.P)))

        self.ts.append(self.t)
        self.t_next_predicts.append(self.t_next_predict)

    def update(self, z, t=numpy.nan):
        if t is numpy.nan:
            self.ukf.update(z)
        else:
            # The update is back dated so we find the time at which it was taken
            index = numpy.searchsorted(self.ts, t) - 1
            ts_old = self.ts[index:]

            # Remove now invalid data
            self._Xs = self._Xs[:index]
            self._Ps = self._Ps[:index]
            self._deviations = self._deviations[:index]
            self.ts = self.ts[:index]
            self.t_next_predicts = self.t_next_predicts[:index]

            # Reset the UKF for sigma calc
            self.ukf.x = self._Xs[-1]
            self.ukf.P = self._Ps[-1]
            self.t = self.ts[-1]
            self.t_next_predict = self.t_next_predicts[-1]
            self.step(ts_old[0] - self.t)

            # Do the update
            self.ukf.update(z)

            # Step forward in time again
            for t_i in ts_old[1:]:
                self.step(t_i - self.t)

    def get_Xs(self):
        return numpy.array(self._Xs)

    def get_deviations(self):
        return numpy.array(self._deviations)

    def get_data(self):
        return numpy.concatenate([self.get_Xs(), self.get_deviations()], axis=1)
