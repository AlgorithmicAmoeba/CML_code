# Contains code for the system model
import numpy


class Model:
    def __init__(self):
        pass

    def DEs(self, X, t, inputs):
        Ng, Nx, Nfa, Ne, Nco, No, Nn, Nb, Nz, Ny, V, Vg = [max(0, N) for N in X]
        Fg_in, Cg_in, Fco_in, Cco_in, Fo_in, Co_in, \
            Fg_out, Cn_in, Fn_in, Fb_in, Cb_in, Fm_in, Fout, Tamb, Q = inputs(t)

        alpha, PO, gamma, theta, beta = 0.1, 1.5, 1.8, 0.0001, 0.1
        delta = 0.2

        # Concentrations
        Cg, Cx, Cfa, Ce, Cn, Cb, Cz, Cy = [N/V for N in [Ng, Nx, Nfa, Ne, Nn, Nb, Nz, Ny]]
        Cco, Co = [N/Vg for N in [Nco, No]]

        # Rate equations:
        # 1) glucose + CO2 + 3*ATP --> 2*FA + 2*water
        # 2) glucose --> 6*CO2 + 12*NADH + 4*ATP (TCA)
        # 3) NADH + 0.5*O2 -> 7/3 ATP (Respiration)
        # 4) glucose -> 2*ethanol + 2*CO2 + 2*ATP
        # 5) glucose + gamma*ATP --> 6*biomass + beta*NADH
        # Unkowns: rFAp, rTCA, rResp, rEp, rXp (5)

        rate_matrix = numpy.array([[1, 0, 0, 0, 0],
                                   [0, 0, 0, 1, 0],
                                   [0, 0, 0, 0, 1],
                                   [-3, 4, 7/3, 2, gamma],
                                   [0, 12, -1, 0, beta]])
        rZ = (0.6 / 46 / 40 * 3 + 2 / 46 / 120 * 3) if Cz > 0 else 0  # decrease
        rY = (0.6 / 46 / 25 * 3 + 0.6 / 46 / 40)*Cy*4  # increase

        rFA_calc = 0.59 * (Cg / (7 + Cg))
        rE_calc = (2/46/120*3.2 + rY - rZ) * (Cg / (1e-3 + Cg))
        RHS = [rFA_calc, rE_calc, 0, theta, 0]

        rFAp, rTCA, rResp, rEp, rXp = numpy.linalg.inv(rate_matrix) @ RHS

        rG = -rFAp - rTCA - rEp - rXp
        rX = 6*rXp
        rFA = 2*rFAp + 0.5*rZ
        rE = 2*rEp
        rCO = -2*rFAp + 6*rTCA + 2*rEp + alpha*rXp
        rO = -0.5*rResp

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
        dNb = Fb_in*Cb_in - Fout*Cb
        dNz = -190*rZ*Cx*V
        dNy = -95*rY*Cx*V
        dV = Fg_in + Fn_in + Fb_in + Fm_in - Fout
        dVg = Fco_in + Fo_in - Fg_out

        return dNg, dNx, dNfa, dNe, dNco, dNo, dNn, dNb, dNz, dNy, dV, dVg

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
        pass

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

