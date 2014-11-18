#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

class mc_packet(object):
    def __init__(self, tau_low):
        self.tau_low = tau_low
        self.tau = tau_low
        self.mu = np.sqrt(np.random.rand(1)[0])

        self.draw_new_tau()

    def draw_new_tau(self):

        self.tau_int = -np.log(np.random.rand(1)[0])

    def scatter_packet(self):

        self.mu = 2. * np.random.rand(1)[0] - 1.

    def advance_packet(self):

        self.tau -= self.mu * self.tau_int

    def check_for_backscatter(self):

        if self.tau > self.tau_low:
            return True
        else:
            return False

    def check_for_escape(self):

        if self.tau < 0:
            return True
        else:
            return False

    def propagate_packet_until_escape(self):

        while True:
            self.advance_packet()
            if self.check_for_backscatter():
                return False
            if self.check_for_escape():
                return True
            self.scatter_packet()
            self.draw_new_tau()


def run_simulation(Npackets, tau_min = 10):

    mu_esc = []
    mu_bins = np.linspace(0, 1, 21)
    dmu = mu_bins[1] - mu_bins[0]

    for i in xrange(Npackets):

        packet = mc_packet(tau_min)
        if packet.propagate_packet_until_escape():

            mu_esc.append(packet.mu)


    mu_esc = np.array(mu_esc)
    mu_mid = (np.floor(mu_esc / dmu) + 0.5 ) * dmu



    plt.hist(mu_esc, bins = mu_bins,  weights = 1. / mu_mid, normed = True, label = "Monte Carlo results")
    plt.plot(mu_bins, 1. / 0.7 * (0.4 + 0.6 * mu_bins), label = "expectation", color = "red", lw = 2)
    plt.xlabel(r"$\mu$")
    plt.ylabel(r"$I(\mu)$ [arbitrary units]")
    plt.legend(loc = "upper left")
    plt.show()


if __name__ == "__main__":

    run_simulation(100000)

