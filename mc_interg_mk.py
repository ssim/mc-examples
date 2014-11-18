import sys
import numpy
from numpy import *

def f(x,n):
    return n/sqrt(pi)*exp(-x**2*n**2)

def int_by_rejection(n,maxiter):
    """ integrattion-by-rejection of f(x,n)
    "hit-or-miss method"
    """

    px=[]
    py=[]

    hits = 0
    trys = 0

    for iter in range(1,maxiter+1):

        # sample from uniform distribution [-1,1]
        x=rand()*2-1
        # sample from uniform distribution [0,n/sqrt(pi)]
        y=rand()*n/sqrt(pi)

        if (y <= f(x,n)):
            hits = hits + 1

        trys = trys + 1

        if (iter % 10000 == 0):

            # ratio of hits to trys
            rat = float(hits)/float(trys)

            # estimated value for the integral
            # (multiply rat by the area of the sampling interval)
            est = rat*2*n/sqrt(pi)

            # estimated variance of rat
            var_rat = rat*(1.-rat)/iter

            # estimated error of est
            err_est = sqrt(var_rat*2*n/sqrt(pi))

            # real error
            error = est-1.

            print est, err_est, error
            sys.stdout.flush()

            px.append(iter)
            py.append(est)

    plot(px,py)

    return float(hits)/float(trys)*2*n/sqrt(pi)


def int_by_mean(n,maxiter):
    """ integrattion-by-mean-value of f(x,n)
    """

    a = -1.0
    b = 1.0

    # define uniform distribution in [a,b]
    def u(x):
        return 1./(b-a)

    px=[]
    py=[]

    meansum = 0.0

    for iter in range(1,maxiter+1):

        # sample from uniform distribution u(x)
        x = (b-a)*rand()+a

        meansum = meansum + f(x,n)/u(x)

        if (iter % 10000 == 0):
            # current estimate of the integral
            est = meansum/iter

            error = est-1.

            print est, error
            sys.stdout.flush()

            px.append(iter)
            py.append(est)


    plot(px,py)

    return meansum/maxiter

def int_by_mean_normal(n,maxiter):
    """ integrattion-by-mean-value of f(x,n)
    with Gaussian importance sampling
    """

    s = 1.0/n

    # define Gaussian distribution consistent with the
    # numpy.random.normal function
    def g(x):
        return exp(-0.5*x**2/s**2)/sqrt(2.*pi*s**2)

    px=[]
    py=[]

    meansum = 0.0

    for iter in range(1,maxiter+1):

        # sample from the distribution function g(x)
        x = numpy.random.normal(0.0,s)

        meansum = meansum + f(x,n)/g(x)

        if (iter % 10000 == 0):

            # current estimate for the integral
            est = meansum/iter

            # absolute error
            error = est-1.

            print est, error
            sys.stdout.flush()

            px.append(iter)
            py.append(est)


    plot(px,py)

    return meansum/maxiter

#example usage:
#
#niter=1000000
#int_by_rejection(5,niter)
#int_by_mean(5,niter)
#int_by_mean_normal(5,niter)







