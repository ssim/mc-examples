#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sp


class homogeneous_plane(object):

    def __init__(self, Lmax, kappa, Npackets):
        self.Lmax = Lmax
        self.kappa = kappa
        self.tau_max = self.Lmax * self.kappa
        self.Npackets = Npackets

        print "Plane Total optical depth: %f" % self.tau_max

    def initialise_packets(self):

        self.xi = self.Lmax * np.random.rand(self.Npackets)
        self.mui = 2. * np.random.rand(self.Npackets) - 1.
        self.tau = -np.log(np.random.rand(self.Npackets))

    def calculate_distance_to_edge(self):

        self.d = np.where(self.mui < 0, -self.xi / self.mui, (self.Lmax - self.xi) / self.mui)
        self.tau_edge = self.d * self.kappa

    def calculate_number_of_absorbed_packets(self):

        self.Nabs = np.sum(self.tau < self.tau_edge)

    def calculate_escape_probability(self):

        self.initialise_packets()
        self.calculate_distance_to_edge()
        self.calculate_number_of_absorbed_packets()

        self.Pesc = 1. - self.Nabs / float(self.Npackets)
        self.Pesc_error = 1. / np.sqrt(float(self.Npackets))

class homogeneous_sphere(object):
    def __init__(self, Rmax, kappa, Npackets):
        self.Rmax = Rmax
        self.kappa = kappa
        self.tau_max = self.Rmax * self.kappa
        self.Npackets = Npackets

        print "Sphere total optical depth: %f" % self.tau_max

    def initialise_packets(self):

        self.ri = self.Rmax * np.random.rand(self.Npackets)**(1./3.)
        self.mui = 2. * np.random.rand(self.Npackets) - 1.
        self.tau = -np.log(np.random.rand(self.Npackets))


    def calculate_distance_to_edge(self):

        self.d = -self.ri * self.mui + np.sqrt(self.ri**2 * (self.mui**2 - 1.) + self.Rmax**2)
        self.tau_edge = self.d * self.kappa

    def calculate_number_of_absorbed_packets(self):

        self.Nabs = np.sum(self.tau < self.tau_edge)

    def calculate_escape_probability(self):

        self.initialise_packets()
        self.calculate_distance_to_edge()
        self.calculate_number_of_absorbed_packets()

        self.Pesc = 1. - self.Nabs / float(self.Npackets)
        self.Pesc_error = 1. / np.sqrt(float(self.Npackets))


def analytic_solution(tau):

    #see Osterbrock - Astrophysics of Gaseous Nebulae, Appendix B
    return 0.75 / tau * (1. - 0.5 / tau**2 + (1. / tau + 0.5 / tau**2) * np.exp(-2. * tau))

def analytic_solution_plane(tau):

    return 0.5 / tau * (1. + np.exp(-tau) * (tau - 1. + tau**2 * np.exp(tau) * sp.expi(-tau)))


if __name__ == "__main__":

    Npackets = 10000

    tau_trial = np.logspace(-3, 2, 1000)
    tau_test = np.array([0.01, 0.1, 0.5, 1, 2, 3, 4, 5, 10, 100])

    Rmax  = 1.
    kappa = tau_test / Rmax

    Pesc = []
    PescE = []

    Pescp = []
    PescEp = []

    for k in kappa:
        sphere = homogeneous_sphere(Rmax, k, Npackets)
        sphere.calculate_escape_probability()

        plane = homogeneous_plane(Rmax, k, Npackets)
        plane.calculate_escape_probability()

        Pesc.append(sphere.Pesc)
        PescE.append(sphere.Pesc_error)

        Pescp.append(plane.Pesc)
        PescEp.append(plane.Pesc_error)

    plt.figure()
    plt.title("homogeneous sphere, %d packets" % sphere.Npackets)
    plt.plot(tau_trial, analytic_solution(tau_trial), label = "analytic solution")
    plt.plot(tau_trial, 0.75 / tau_trial, color = "blue", ls = "dashed", label = r"$3/4\tau$")
    plt.errorbar(tau_test, Pesc, yerr = PescE, ls = "", color = "red", marker = "o", label = "Monte Carlo results")
    plt.ylim([1e-2,1.1])
    plt.xlabel(r"$\tau$")
    plt.ylabel(r"$P_{\mathrm{esc}}$")
    plt.legend(loc = "lower left")
    plt.xscale("log")
    plt.yscale("log")

    plt.figure()
    plt.title("homogeneous slab, %d packets" % sphere.Npackets)
    plt.plot(tau_trial, analytic_solution_plane(tau_trial), color = "blue", label = "analytic solution")
    plt.plot(tau_trial, 0.5 / tau_trial, color = "blue", ls = "dashed", label = r"$1/2\tau$")
    plt.errorbar(tau_test, Pescp, yerr = PescEp, ls = "", color = "red", marker = "o", label = "Monte Carlo results")
    plt.ylim([1e-2,1.1])
    plt.xlabel(r"$\tau$")
    plt.ylabel(r"$P_{\mathrm{esc}}$")
    plt.xscale("log")
    plt.yscale("log")

    plt.legend(loc = "lower left")

    plt.show()

