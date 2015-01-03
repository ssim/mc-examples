#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt

c = 2.9979e10

class mc_packet(object):
    def __init__(self, Rmin, Rmax, nu_min, nu_max, nu_line, tau_sob, t, verbose = False):
        self.verbose = verbose
        self.Rmin = Rmin
        self.Rmax = Rmax

        self.nu_line = nu_line
        self.tau_sob = tau_sob

        self.r = self.Rmin
        self.mu = np.sqrt(np.random.rand(1)[0])
        self.nu = nu_min + (nu_max - nu_min) * np.random.rand(1)[0]

        self.t = t
        self.emergent_nu = []

        self.draw_new_tau()
        self.check_for_boundary_intersection()
        self.calc_distance_to_sobolev_point()

    def draw_new_tau(self):

        self.tau_int = -np.log(np.random.rand(1)[0])

    def update_position_direction(self, l):

        ri = self.r

        self.r = np.sqrt(self.r**2 + l**2 + 2 * l * self.r * self.mu)
        self.mu = (l + self.mu * ri) / self.r


    def check_for_boundary_intersection(self):

        if self.mu <= -np.sqrt(1 - (self.Rmin / self.r)**2):
            ## packet will intersect inner boundary if not interrupted
            sgn = -1.
            rbound = self.Rmin
            self.boundint = "left"
        else:
            ## packet will intersect outer boundary if not interrupted
            sgn = 1.
            rbound = self.Rmax
            self.boundint = "right"

        self.lbound = -self.mu * self.r + sgn * np.sqrt((self.mu * self.r)**2 - self.r**2 + rbound**2)


    def perform_interaction(self):

        mui = self.mu
        beta = self.r / self.t / c

        self.mu = 2. * np.random.rand(1)[0] - 1.
        self.mu = (self.mu + beta) / (1 + beta * self.mu)

        self.nu = self.nu * (1. - beta * mui) / (1. - beta * self.mu)


    def calc_distance_to_sobolev_point(self):

        self.lsob = c * self.t * (1 - self.nu_line / self.nu) - self.r * self.mu

    def propagate(self):

        while True:
            if self.lbound < self.lsob:
                if self.verbose:
                    print("Reaching boundary")
                if self.boundint == "left":
                    if self.verbose:
                        print("Intersecting inner boundary")
                    self.emergent_nu = None
                    break
                else:
                    if self.verbose:
                        print("Escaping through outer boundary")
                    self.emergent_nu = self.nu
                    break
            else:
                if self.verbose:
                    print("Reaching Sobolev point")
                self.update_position_direction(self.lsob)
                if self.tau_sob >= self.tau_int:
                    if self.verbose:
                        print("Line Interaction")
                    self.perform_interaction()


            self.draw_new_tau()
            self.check_for_boundary_intersection()
            self.calc_distance_to_sobolev_point()


class homologous_sphere(object):
    def __init__(self, Rmin, Rmax, nu_min, nu_max, nu_line, tau_sob, t, npack,  verbose = False):

        self.Rmin = Rmin
        self.Rmax = Rmax

        self.nu_min = nu_min
        self.nu_max = nu_max

        self.nu_line = nu_line
        self.tau_sob = tau_sob

        self.t = t

        self.verbose = verbose
        self.packets = [mc_packet(Rmin, Rmax, nu_min, nu_max, nu_line, tau_sob, t, verbose = verbose) for i in xrange(npack)]

        self.emergent_nu = []

    def perform_simulation(self):

        for pack in self.packets:
            pack.propagate()
            if pack.emergent_nu is not None:
                self.emergent_nu.append(pack.emergent_nu)

        self.emergent_nu = np.array(self.emergent_nu)


def main():

    lam_line = 1215.6 * 1e-8
    lam_min = 1000.0 * 1e-8
    lam_max = 1400.0 * 1e-8
    tau_sob = 0.1

    t = 13.5 * 86400.

    vmin = 1e-4 * c
    vmax = 0.01 * c

    Rmin = vmin * t
    Rmax = vmax * t

    nu_min = c / lam_max
    nu_max = c / lam_min
    nu_line = c / lam_line

    npack = 10
    verbose = True

    sphere = homologous_sphere(Rmin, Rmax, nu_min, nu_max, nu_line, tau_sob, t, npack,  verbose = verbose)
    sphere.perform_simulation()

    print(sphere.emergent_nu)

if __name__ == "__main__":

    main()



