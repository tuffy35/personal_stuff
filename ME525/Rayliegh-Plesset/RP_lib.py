import numpy as np
from numpy import pi, sin, array
import matplotlib.pyplot as plt
import pickle

def sphere_volume(r):
    return (4./3.) * pi * r**3

class Bubble(object):
    def __init__(self):
        self.R0 = 0.5e-4
        self.R_start = 0.5e-4
        self.dR_start = 0.
        self.p0 = 101325.
        self.pvap = 0.093 * 1e5
        self.s = (25.2 + 25 * 0.1120) * 1e-3
        self.rho = 790.
        self.nu = 0.4 * 1e-3
        self.c = 1174.
        self.kappa = 1.4
        self.amp = 0.6 * self.p0
        self.f = 2e4
        self.t = []
        self.R = []
        self.Pg = []
        self.t_start = 0.
        self.dt_coarse = 1e-9
        self.dt_fine = 1e-12
        self.t_run = 19e-6
        self.bubble_radiates = True

    def p_acoustic(self, t):
        return -self.amp * sin(2 * pi * self.f * t)

    def dR(self, r, t):
        R, dR = r
        p_gas = (self.p0 + 2*self.s/self.R0 - self.pvap) * (self.R0/R)**(3*self.kappa)
        p_surf = 2 * self.s / R
        p_liq = p_gas + self.pvap - p_surf
        p_ext = self.p0 + self.p_acoustic(t)
        if self.bubble_radiates and len(self.t) > 2:
            Pgdot = (self.Pg[-1] - self.Pg[-2]) / (self.t[-1] - self.t[-2])
            radiation_loss = Pgdot * R / self.c
            ddR = -3*dR**2/(2*R) + 1/(self.rho*R) * (p_liq - 4*self.nu*dR/R - p_ext + radiation_loss)
        else:
            ddR = -3*dR**2/(2*R) + 1/(self.rho*R) * (p_liq - 4*self.nu*dR/R - p_ext)
        return np.array([dR, ddR])

    def calculate_Pgas(self, r):
        return (self.p0 + 2*self.s/self.R0 - self.pvap) * (self.R0/r)**(3*self.kappa)

    def integrate_RK4(self):
        t = self.t
        w = self.R
        t.append(self.t_start)
        w.append(array([self.R_start, self.dR_start]))
        self.Pg.append(self.calculate_Pgas(self.R_start))
        while t[-1] < self.t_start + self.t_run:
            if w[-1][0] < 0.6 * self.R0:
                dt = 0.02 * self.dt_fine if w[-1][0] < 0.06 * self.R0 else self.dt_fine
            else:
                dt = self.dt_coarse
            t.append(t[-1] + dt)
            h = t[-1] - t[-2]
            k1 = self.dR(w[-1], t[-2])
            k2 = self.dR(w[-1] + 0.5*h*k1, t[-2] + 0.5*h)
            k3 = self.dR(w[-1] + 0.5*h*k2, t[-2] + 0.5*h)
            k4 = self.dR(w[-1] + h*k3, t[-1])
            w.append(w[-1] + (h*k1 + 2*h*k2 + 2*h*k3 + h*k4)/6.)
            self.Pg.append(self.calculate_Pgas(w[-1][0]))

    def integrate_RK3(self):
        t = self.t
        w = self.R
        t.append(self.t_start)
        w.append(array([self.R_start, self.dR_start]))
        self.Pg.append(self.calculate_Pgas(self.R_start))

        while t[-1] < self.t_start + self.t_run:
            if w[-1][0] < 0.6 * self.R0:
                dt = 0.02 * self.dt_fine if w[-1][0] < 0.06 * self.R0 else self.dt_fine
            else:
                dt = self.dt_coarse

            t.append(t[-1] + dt)
            h = t[-1] - t[-2]

            k1 = self.dR(w[-1], t[-2])
            k2 = self.dR(w[-1] + 0.5* h*k1, t[-2] + 0.5*h)
            k3 = self.dR(w[-1] - h*k1 + 2*h*k2, t[-2] + h)

            w.append(w[-1] + (h / 6.0) * (k1 + 4 * k2 + k3))
            self.Pg.append(self.calculate_Pgas(w[-1][0]))


    def list2array(self):
        self.t = array(self.t)
        self.R = array(self.R)
        self.Pg = array(self.Pg)

    def pickle_self(self, suffix):
        with open('pickled_bubble' + suffix + '.txt', 'wb') as f:
            pickle.dump(self, f)

def unpickle_thing(path):
    with open(path, 'rb') as f:
        return pickle.load(f)

