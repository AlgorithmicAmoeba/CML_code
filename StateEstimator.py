import numpy
import scipy
import filterpy.kalman
from AdjMerweScaledSigmaPoints import MerweScaledSigmaPoints
import model


class StateEstimator:
    def __init__(self, X0, inputs, t_predict):
        self.inputs = inputs

        self._Xs = [X0]
        self._Ps = [[0]*len(X0)]

        #                           Ng, Nx, Nfa, Ne, Nco, No, Nn, Na, Nb, Nz, Ny, V, Vg
        self.Q = numpy.diag(numpy.array([1e-4, 1e-3, 1e-5, 1e-1, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-2, 1e-2, 1e-5, 1e-5]))
        self.R = numpy.diag(numpy.array([1e-12, 1e-12, 1e-12]))

        self.fx = self.FXObj(self.inputs)
        self.nx = len(self.Q)

        self.sigmas = MerweScaledSigmaPoints(self.nx, 1e-3, 2, 0, sqrt_method=scipy.linalg.sqrtm)
        self.ukf = filterpy.kalman.UnscentedKalmanFilter(self.nx, 3, 0, self.hx, self.fx, self.sigmas)

        self.ukf.x = X0
        self.ukf.Q = self.Q
        self.ukf.R = self.R

        self.t_next_predict = 0
        self.t_predict = t_predict

    @staticmethod
    def hx(x):
        Ng, _, Nfa, Ne, _, _, _, _, _, _, _, V, _ = x
        return Ng/V, Nfa/V, Ne/V

    class FXObj:
        def __init__(self, inputs, t=0):
            self.t = t
            self.inputs = inputs

        def __call__(self, x, dt):
            ts = numpy.linspace(0, dt, 100)
            dt_small = ts[1]
            m_f = model.Model(x, self.inputs)
            m_f.t = self.t
            for _ in ts:
                m_f.step(dt_small)
            return m_f.X

    def step(self, t):
        if t > self.t_next_predict:
            self.ukf.predict(self.t_predict)
            self.t_next_predict += self.t_predict

        self._Xs.append(self.ukf.x)
        self._Ps.append(numpy.sqrt(numpy.diag(self.ukf.P)))

    def update(self, z, t):
        self.ukf.update(z)
        self.fx.t = t

    def get_Xs(self):
        return numpy.array(self._Xs)

    def get_Ps(self):
        return numpy.array(self._Ps)
