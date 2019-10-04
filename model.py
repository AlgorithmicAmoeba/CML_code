# Contains code for the system model
import numpy


class Model:
    def __init__(self, X0):
        self.X = X0
        self.t = 0

    def DEs(self, t, inputs):
        """
        Contains the differential and algebraic equations for the system model.
        The rate equations defined in the matrix `rate_matrix` are described by:
        1) glucose + 2*CO2 + 6*ATP --> 2*FA + 2*water
        2) glucose --> 6*CO2 + 12*NADH + 4*ATP (TCA)
        3) NADH + 0.5*O2 -> 7/3 ATP (Respiration)
        4) glucose -> 2*ethanol + 2*CO2 + 2*ATP
        5) glucose + gamma*ATP --> 6*biomass + beta*NADH
        where the unknowns are: rFAp, rTCA, rResp, rEp, rXp

        Parameters
        ----------
        t : int
            The current Time

        inputs : array_like
            The inputs to the system at the current time

        Returns
        -------
        dX : array_like
            The differential changes to the state variables
        """
        Ng, Nx, Nfa, Ne, Nco, No, Nn, Na, Nb, Nz, Ny, V, Vg = [max(0, N) for N in self.X]
        Fg_in, Cg_in, Fco_in, Cco_in, Fo_in, Co_in, \
            Fg_out, Cn_in, Fn_in, Fb_in, Cb_in, Fm_in, Fout, Tamb, Q = inputs(t)

        alpha, PO, gamma, theta, beta = 0.1, 0.1, 1.8, 0.1, 0.1
        delta = 0.2

        # Concentrations
        Cg, Cx, Cfa, Ce, Cn, Ca, Cb, Cz, Cy = [N/V for N in [Ng, Nx, Nfa, Ne, Nn, Na, Nb, Nz, Ny]]
        Cco, Co = [N/Vg for N in [Nco, No]]

        rate_matrix = numpy.array([[1, 0, 0, 0, 0],
                                   [0, 0, 0, 1, 0],
                                   [0, 0, 0, 0, 1],
                                   [-6, 4, 7/3, 2, -gamma],
                                   [0, 12, -1, 0, beta]])
        first_increase = (0.6 / 46 / 25 * 4) * Cy * 1.8
        second_increase = 2/46/120*3.2
        decrease = (0.6 / 46 / 40*3) * Cz/3
        rZ = decrease + second_increase if Cz > 0 else 0  # decrease
        rY = first_increase + decrease if Cy > 0 else 0  # increase

        rFA_calc = 15e-3 * (Cg / (1e-2 + Cg)) - 0.5*rZ
        rE_calc = (second_increase + rY) * (Cg / (1e-5 + Cg))
        theta_calc = theta * (Cg / (1e-3 + Cg))
        RHS = [rFA_calc, rE_calc, 8e-5, theta_calc, 0]

        rFAp, rTCA, rResp, rEp, rXp = numpy.linalg.inv(rate_matrix) @ RHS

        rG = -rFAp - rTCA - rEp - rXp
        rX = 6*rXp
        rFA = 2*(rFAp + 0.5*rZ)
        rE = 2*(rEp - rZ) * (Cg / (1e-5 + Cg))
        rCO = -2*rFAp + 6*rTCA + 2*rEp + alpha*rXp
        rO = -0.5*rResp
        print(rG)

        # pH calculations
        # Kna, Kfa = 10 ** (14 - 0.2), 10 ** 3.03
        # Ch = Cfa_m - Coh
        # Cfa = Cfa_m + Cfa_u
        # Cb = Cna + Coh
        # Kfa = Ch * Cfa_m / Cfa
        # Kna = Cb * Ch / Cna
        # solving in sympy gives:
        # Ch = (Cb*Kna - numpy.sqrt(Kna*(Cb**2*Kna - 4*Cb*Cfa*Kfa + 4*Cfa*Kfa*Kna)))/(2*(Cb - Kna))
        # pH = -numpy.log10(Ch)

        # DE's
        dNg = Fg_in*Cg_in - Fout*Cg + rG*Cx*V
        dNx = rX*Cx*V
        dNfa = -Fout*Cfa + rFA*Cx*V
        dNe = -Fout*Ce + rE*Cx*V
        dNco = Fco_in*Cco_in - Fg_out*Cco + rCO*Cx*V
        dNo = Fo_in*Co_in - Fg_out*Co - rO*Cx*V
        dNn = Fn_in*Cn_in - Fout*Cn - delta*rX*Cx*V
        dNa = - Fout * Ca
        dNb = Fb_in*Cb_in - Fout*Cb
        dNz = -190*rZ*Cx*V
        dNy = -95*rY*Cx*V
        dV = Fg_in + Fn_in + Fb_in + Fm_in - Fout
        dVg = Fco_in + Fo_in - Fg_out

        return dNg, dNx, dNfa, dNe, dNco, dNo, dNn, dNa, dNb, dNz, dNy, dV, dVg

    def step(self, inputs, dt):
        """
        Updates the model with inputs
        Parameters
        ----------
        inputs : array_like
            List of all the inputs to the model
        dt : float
            Amount of time since the previous step

        Returns
        -------
        outputs : array_like
            List of all the outputs from the model
        """
        self.t += dt
        dX = self.DEs(self.t, inputs)
        self.X += numpy.array(dX)*dt
        return self.X

    def instrument_bias_gain(self, raw_inputs):
        """
        Takes in the raw inputs and translates them into their
        actual values using the gain and bias
        Parameters
        ----------
        raw_inputs : array_like
            List of raw input values

        Returns
        -------
        inputs : array_like
            List of actual input values

        """
        pass

