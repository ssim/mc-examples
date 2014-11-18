#!/usr/bin/env python
import numpy as np
import numpy.random as rng
import scipy.integrate as integ
import matplotlib.pyplot as plt


def f(x, n):

    return n / np.sqrt(np.pi) * np.exp(-(n * x)**2)


def find_appropriate_n(nmax = 10, acc = 1e-5):

    ntrial = np.arange(1, nmax)

    for ni in ntrial:

        integral = integ.quad(f, -1, 1, args = (ni,))[0]

        if np.fabs(integral - 1.) < 1e-5:

            print "integral for n = %d within requested accuracy of %e" % (ni, acc)
            nsol = ni
            break
    else:
        print "Warning: Tested until %d: No appropriate n found" % nmax
        nsol = -1

    return nsol

def integration_by_rejection(Nsamples, nsol):

    ymax = f(0, nsol)

    points = rng.rand(2, Nsamples)

    x = 2. * points[0] - 1.
    y = ymax * points[-1]

    nbelow = np.sum(y <= f(x, nsol))

    integral = nbelow / float(Nsamples) * 2. * ymax

    return integral

def integration_by_mean_value(Nsamples, nsol):

    x = rng.rand(Nsamples)

    fmean = np.mean(f(x, nsol))

    integral = 2. * fmean

    return integral

def integration_by_mean_value_trans(Nsamples, nsol):

    g = lambda x, sigma, lam: 1. / (sigma * np.sqrt(2. * np.pi)) * np.exp(-0.5 * ((x - lam) / sigma)**2)
    x = rng.standard_normal(Nsamples)

    sigma = 1.
    lam = 0.

    fmean = np.mean(f(x, nsol) / g(x, sigma, lam))

    integral = fmean

    return fmean






if __name__ == "__main__":

    nsol = find_appropriate_n()
    assert(nsol != -1)

    print integration_by_rejection(1000., nsol)
    print integration_by_mean_value(1000., nsol)
    print integration_by_mean_value_trans(1000., nsol)

