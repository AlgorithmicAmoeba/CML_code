# Contains code for the system model
import numpy


class Model:
    def __init__(self):
        pass

    def DEs(self, X, t, inputs):
        Ng, Nx, Nfa, Ne, Nco, No, Nn, Nb, V, Vg = [max(0, N) for N in X]
        Fg_in, Cg_in, Fco_in, Cco_in, Fo_in, Co_in, \
            Fg_out, Cn_in, Fn_in, Fb_in, Cb_in, Fm_in, Fout, Tamb, Q = inputs(t)
        a, b, c = 1, 0, 10
        # a: how much does Cg affect rG
        # b: how much does Cn affect rX
        # c: ratio of rFA to rE: rFA = c*rE
        alpha, PO, gamma, theta, beta = 0.1, 1.5, 1.8, 0.0001, 0.1
        delta = 0.2

        # Concentrations
        Cg, Cx, Cfa, Ce, Cn, Cb = [N/V for N in [Ng, Nx, Nfa, Ne, Nn, Nb]]
        Cco, Co = [N/Vg for N in [Nco, No]]

        # Rate equations:
        # 1) rG = a * Cg
        # 2) rX = b * Cn
        # 3) rFA = c * rE
        # 4) rG = rPyr + (1 + alpha)*rX
        # 5) rPyr = 0.75*rFA  + rCO + 1.5*rE
        # 6) 1/3*rPyr + 1/3*rCO + 2*PO*rO - gamma*rX - 0.25*rFA = theta : ATP
        # 7) 0.25*rFA + 0.5*rE + 2*rO - beta*rX - 1/3*rPyr - 5/3*rCO = 0 : NADH
        # Unkowns: rG, rX, rFA, rE, rPyr, rCO, rO (7)

        rate_matrix = numpy.array([[0, 0, 1, 0, 0, 0, 0],
                                   [0, 1, 0, 0, 0, 0, 0],
                                   [0, 0, 0, 0, 0, 1, 0],
                                   [-1, 1 + alpha, 0, 0, 1, 0, 0],
                                   [0, 0, 0.75, 1.5, -1, 1, 0],
                                   [0, -gamma, -0.25, 0, 1/3, 1/3, 2*PO],
                                   [0, beta, -0.25, -0.5, 1/3, 5/3, -2]])
        rFA_calc = 0.001429825802986 * 1.35  # (Cg / (1e-3 + Cg))
        RHS = [rFA_calc, b*Cn, 0, 0, 0, theta, 0]

        rG, rX, rFA, rE, rPyr, rCO, rO = numpy.linalg.inv(rate_matrix) @ RHS
        print(rG, rX, rFA, rE, rPyr, rCO, rO)
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
        dNg = Fg_in*Cg_in - Fout*Cg - 6*rG*Cx*V
        dNx = rX*Cx*V
        dNfa = -Fout*Cfa + 4*rFA*Cx*V
        dNe = -Fout*Ce + 2*rE*Cx*V
        dNco = Fco_in*Cco_in - Fg_out*Cco + rCO*Cx*V
        dNo = Fo_in*Co_in - Fg_out*Co - rO*Cx*V
        dNn = Fn_in*Cn_in - Fout*Cn - delta*rX*Cx*V
        dNb = Fb_in*Cb_in - Fout*Cb
        dV = Fg_in + Fn_in + Fb_in + Fm_in - Fout
        dVg = Fco_in + Fo_in - Fg_out

        return dNg, dNx, dNfa, dNe, dNco, dNo, dNn, dNb, dV, dVg

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

