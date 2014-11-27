#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt


c = 2.9979e10

class packet(object):
    def __init__(self, chi = 1):
        self.x = 0
        self.mu = 2. * np.random.rand(1)[0] - 1.
        self.chi = chi
        self.lacc = 0

        self.calculate_interaction_location()

    def calculate_interaction_location(self):
        self.lint = -np.log(np.random.rand(1)[0]) / self.chi

    def propagate_packet(self, l):

        self.x = self.x + self.mu * l
        self.lacc = self.lacc + l

    def scatter_packet(self):

        self.mu = 2. * np.random.rand(1)[0] - 1.

    def propagate(self, tend):

        lend = tend * c
        while (self.lacc + self.lint) < lend:

            self.propagate_packet(self.lint)
            self.scatter_packet()
            self.calculate_interaction_location()

        self.propagate_packet(lend - self.lacc)
        self.lint -= (lend - self.lacc)


class packet_pool(object):
    def __init__(self, N, chi = 1):
        self.N = int(N)
        self.packets = [packet(chi = chi) for i in xrange(self.N)]
        self.t = 0

    def propagate_step(self, t):

        [self.packets[i].propagate(t) for i in xrange(self.N)]
        self.t += t

    def return_locations(self):

        return np.array([self.packets[i].x for i in xrange(self.N)]).copy()

    def propagate_steps(self, tsteps):

        xs = {}

        for t in tsteps:
            print("propagate until ", t)
            self.propagate_step(t)
            xs[t] = self.return_locations()

        return xs


def analytic_solution(x, t,  L, Etot, chi):

    D = c / 3. / chi

    return Etot / np.sqrt(4. * np.pi * D * t) * np.exp(-x**2 / (4 * D *t))


if __name__ == "__main__":

    np.random.seed(12314)
    Lend = 1
    tend = 0.5 * Lend / c * 10
    tend = 5e-12
    Npackets = 1e4
    Etot = 1e10 / 101.
    chi = 1e2
    bins = np.linspace(-Lend * 0.5 , 0.5 * Lend, 101)
    dbin = bins[1] - bins[0]

    t = [5e-12, 1e-11, 2e-11, 5e-11]

    pool = packet_pool(Npackets, chi = chi)
    xs = pool.propagate_steps(t)
    x = np.linspace(-Lend * 0.5 , Lend * 0.5, 1000)
    labels = [r"$t = %.1f\times 10^{-11}\,\mathrm{s}$" % (ti * 1e11) for ti in t]

    for i, k in enumerate(sorted(xs.keys())):
        ret = plt.hist(xs[k], bins = bins, histtype = "step", weights = np.ones(Npackets) * Etot / Npackets / dbin, label = labels[i] )

    print(np.sum(ret[0] * (ret[1][1:] - ret[1][:-1])))

    p = [plt.plot(x, analytic_solution(x, ti, Lend, Etot, chi), ls = "solid", color = "black") for ti in t]

    p[0][0].set_label("analytic solution")

    plt.legend()
    plt.xlim([-0.5 * Lend, 0.5 * Lend])
    plt.xlabel(r"$x$")
    plt.ylabel(r"$E$ [$\mathrm{erg\,cm^{-1}}$]")
    plt.savefig("diffusion_results.pdf")

    plt.show()

