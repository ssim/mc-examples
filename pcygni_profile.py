#!/usr/bin/env python
## Important References:
##     Jeffery & Branch 1990: "Analysis of Supernova Spectra"
##     ADS link:http://adsabs.harvard.edu/abs/1990sjws.conf..149J
## 
##     Thomas et al 2011: "SYNAPPS: Data-Driven Analysis for Supernova Spectroscopy"
##     ADS link:http://adsabs.harvard.edu/abs/2011PASP..123..237T
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integ

c = 2.99792458e10

class homologous_sphere(object):
    """
    Class describing a sphere in homologous expansions, i.e. the supernova ejecta
    """
    def __init__(self, rmin = 1e14, rmax = 1e15, vmax = 1e8, Ip = 1, tauref = 1, vref = 5e7, ve = 5e7, lam0 = 1215.7 * 1e-8):
        """
        Keyword arguments:
        rmin -- the inner edge of the ejecta, i.e. the location of the
                photosphere (default 1e14)
        rmax -- the outer edge of the ejecta (default 1e15)
        vmax -- the expansion velocity at the outer edge of the ejecta
                (default 1e8)
        Ip   -- measure for the incident intensity at the photosphere
                (default 1)
        tauref -- Sobolev optical depth at the location vref in the ejecta
                  (default 1)
        vref -- location in the ejecta where the reference optical depth tauref
                is measured (default 5e7)
        ve   -- additional parameter for the density and in turn optical depth
                law in the ejecta, c.f. Thomas et al. 2011 (default 5e7)
        lam0 -- rest frame wavelength (in Angstrom!) of the line transition
                (default 1215.7, i.e. Lyman-alpha)
        """

        self._t = None
        self._nu0 = None
        self._zmax = None

        self.rmin = rmin
        self.rmax = rmax
        self.vmax = vmax
        self.Ip = Ip
        self.tauref = tauref


        self.lam0 = lam0
        self.vref = vref
        self.ve = ve


    @property
    def t(self):
        """Time since explosion; determined from rmax and vmax"""
        if self._t is None:
            self._t = self.rmax / self.vmax
        return self._t

    @property
    def nu0(self):
        """Rest-frame frequency of line transition; determined from lam0"""
        if self._nu0 is None:
            self._nu0 = c / self.lam0
        return self._nu0

    @property
    def zmax(self):
        """Maximum z value in the ejecta (c.f. Jeffrey & Branch 1990); determined from rmax"""
        if self._zmax is None:
            self._zmax = self.rmax
        return self._zmax

    def calc_z(self, nu):
        """
        Calculate location (in terms of z) of resonance point for photon emitted by the photosphere with frequency nu
        
        Arguments:
        nu -- photospheric frequency of photon

        Returns:
        z -- z-location of resonance point
        """

        return c * self.t * (1. - self.nu0 / nu)

    def calc_p(self, r, z):
        """
        Calculate p-coordinate of location (r,z) in ejecta; c.f. Jeffrey &
        Branch 1990 for basic impact geometry in the elementary supernova
        model.

        Arguments:
        r -- radial coordinate of location of interest
        z -- z-coordinate (i.e. along the line-of-sight to the observer) of the location of interest

        Returns:
        p -- p-coordinate (perpendicular to z) of the location of interest
        """

        assert(np.fabs(r) > np.fabs(z))

        return np.sqrt(r**2 - z**2)

    def calc_r(self, p, z):
        """
        Calculate radius of location (z, p) in ejecta;

        Arguments:
        p -- p-coordinate (perpendicular to line-of-sight to observer)
        z -- z-coordinate (along line-of-sight to observer)

        Returns:
        r -- radius of location
        """

        return np.sqrt(p**2 + z**2)

    def calc_W(self, r):
        """
        Calculate geometric dilution factor

        Arguments:
        r -- radius of location

        Returns:
        W -- geometric dilution factor
        """

        return 0.5 * (1. - np.sqrt(1. - (self.rmin / r)**2))

    def calc_tau(self, r):
        """
        Calculate line optical depth at radius r, according to density profile.

        We assume an exponential density and thus optical depth profile as presented
        in Thomas et al. 2011.

        Arguments:
        r -- radius of location

        Returns:
        tau -- line optical depth
        """

        v = r / self.t

        return self.tauref * np.exp((self.vref - v) / self.ve)

    def S(self, p, z, mode = "both"):
        """
        Calculate source function at location (p, z) in ejecta.

        In case only the pure absorption component of the line profile is
        considered, the source function is of course 0. Otherwise, it follows
        from eq. of Jeffery & Branch 1990.

        Arguments:
        p -- p-coordinate of location
        z -- z-coordinate of location
        
        Keyword arguements:
        mode -- whether the source function should be calculated for the pure
                absorption case (mode = 'abs') or including both absorption and
                emission (mod = 'both'); (default = 'both')

        Returns:
        S -- source function at location (p, z)
        """

        if mode == "abs":
            return 0

        r = self.calc_r(p, z)

        if r > self.rmax or r < self.rmin:
            #outside ejecta or inside photosphere"
            return 0
        elif z < 0 and p < self.rmin:
            #occulted region"
            return 0
        else:
            #emission region"
            return self.calc_W(r) * self.Ip

    def I(self, p, z):
        """
        Determine the initial specific intensity at location (p, z)

        Used in eq. of Jeffery & Branch 1990. Only if the line of sight going through (p, z) and towards the observer
        intersects the photosphere, a non-zero initial specific intensity is found.

        Arguments:
        p -- p-coordinate of location of interest
        z -- z-coordinate of location of interest

        Returns:
        I -- initial specific intensity
        """

        if p < self.rmin:
            #in the photosphere plane"
            return self.Ip
        else:
            #above the photosphere plane"
            return 0

    def tau(self, p, z):
        """
        Determine the line optical on the line-of-sight towards the observer,
        at location (p, z).

        Used in eq. of Jeffery & Branch 1990. Only locations in the emission
        region outside of the occulted zone may attenuated the radiation field.
        Thus, only there a non-zero optical depth is returned.

        Arguments:
        p -- p-coordinate of the location of interest
        z -- z-coordinate of the location of interes

        Returns:
        tau -- optical depth at the location of interest
        """

        r = self.calc_r(p, z)

        if r > self.rmax or r < self.rmin:
            #outside ejecta or inside photosphere"
            return 0
        elif z < 0 and p < self.rmin:
            #occulted region"
            return 0
        else:
            #emission region"
            return self.calc_tau(r)


    def Iemit(self, p, z, mode = "both"):
        """
        Determine the total specific intensity at location (p, z).

        The absorption or emission-only cases may be treated, or both effects
        may be included to calculate the full line profile. Used in eq. of Jeffery & Branch 1990.

        Arguments:
        p -- p-coordinate of location of interest
        z -- z-coordinate of location of interest

        Keyword arguments:
        mode -- identifies the included interaction channels: 'abs' for pure
                absorption, 'emit' for pure emission, 'both' for the inclusion
                of both effects (default 'both')

        Returns:
        Itot -- total specific intensity at location (p, z)
        """
        tau = self.tau(p, z)

        if mode == "both" or mode == "abs":
            return (self.I(p, z) * np.exp(-tau) + self.S(p, z, mode = mode) * (1. - np.exp(-tau))) * p
        else:
            return (self.I(p, z) + self.S(p, z) * (1. - np.exp(-tau))) * p

    def calc_line_flux(self, nu, mode = "both"):
        """
        Calculate the emergent line flux at the frequency nu

        Arguments:
        nu -- lab frame frequency at which the line flux is to be calculated

        Keyword arguments:
        mode -- identifies the included interaction channels; see self.Iemit (default 'both')

        Returns:
        Fline -- line flux
        """

        z = self.calc_z(nu)
        pmax = self.rmax

        #integration over impact parameter p
        res = 2. * np.pi * integ.quad(self.Iemit, 0, pmax, args = (z, mode))[0]
        return res

    def calc_line_profile(self, nu_min, nu_max, npoints = 100, mode = "both"):
        """
        Calculate the full line profile between the limits nu_min and nu_max

        Arguments:
        nu_min -- lower frequency limit
        nu_max -- upper frequency limit

        Keyword arguments:
        npoints -- number of points of the equidistant frequency grid (default 100)
        mode    -- identifier setting the interaction mode, see self.Iemit (default 'both')

        Returns:
        nu -- frequency points
        F  -- associated flux
        """

        nu = np.linspace(nu_min, nu_max, npoints)

        F = []

        for nui in nu:
            F.append(self.calc_line_flux(nui, mode = mode))

        return nu, np.array(F)


    def show_line_profile(self, nu_min, nu_max, npoints = 100, include_abs = True, include_emit = True, vs_nu = True):
        """
        Visualise Line Profile

        The P-Cygni line profile will always be displayed. The pure absorption
        and emission components can be included in the plot as well. The flux
        (will always be be F_nu) may be plotted against frequency or
        wavelength.

        Arguments:
        nu_min  -- lower frequency limit
        nu_max  -- upper frequency limit

        Keyword arguments:
        npoints -- number of points of the frequency grid (default 100)
        include_abs  -- if True, the pure absorption flux will be included and
                        shown as a separate line (default True)
        include_emit -- if True, the pure emission flux will be included and
                        shown as a separate line (default True)
        vs_nu -- if True the quantities will be shown against frequency,
                 otherwise against wavelength (default True)

        Returns:
        fig -- figure instance containing plot
        """

        nu, Fline = self.calc_line_profile(nu_min, nu_max, npoints = npoints)

        if include_abs:
            Fabs = self.calc_line_profile(nu_min, nu_max, npoints = npoints, mode = "abs")[-1]

        if include_emit:
            Femit = self.calc_line_profile(nu_min, nu_max, npoints = npoints, mode = "emit")[-1]

        if vs_nu:
            x = nu
        else:
            x = c / nu * 1e8


        fig = plt.figure()
        ax = fig.add_subplot(111)
        fig.subplots_adjust(top = 0.8)

        if include_abs:
            ax.plot(x, Fabs, color = "grey", ls = "dashed", label = "absorption component")
        if include_emit:
            ax.plot(x, Femit, color = "grey", ls = "dotted", label = "emission component")

        ax.plot(x, Fline, color = "blue", ls = "solid", label = "emergent line profile")
        ax.legend(bbox_to_anchor=(0., 1.05, 1., .102), loc=3,ncol=2, mode="expand", borderaxespad=0.)

        if vs_nu:
            ax.set_xlabel(r"$\nu$ [Hz]")
        else:
            ax.set_xlabel(r"$\lambda$ [$\AA$]")

        ax.set_ylabel(r"$F_{\nu}$ [$\mathrm{erg\,cm^{-2}\,s^{-1}\,Hz^{-1}}$]")
        ax.set_xlim([np.min(x), np.max(x)])

        return fig

    def save_line_profile(self, nu_min, nu_max, npoints = 100, include_abs = True, include_emit = True, vs_nu = True, fname = "line_profile.txt"):
        """
        Save Line Profile data to a text file

        The P-Cygni line profile will always be stored. The pure absorption and
        emission components can be included into the output as well. In the
        first column of the data table, either the frequency or the wavelength
        may be stored. 

        Arguments:
        nu_min  -- lower frequency limit
        nu_max  -- upper frequency limit

        Keyword arguments:
        npoints -- number of points of the frequency grid (default 100)
        include_abs  -- if True, the pure absorption flux will be stored in the
                        output as well (default True)
        include_emit -- if True, the pure emission flux will be stored in the
                        output as well (default True)
        vs_nu -- if True the frequency will be stored in the first column,
                 otherwise the wavelength (default True)
        fname -- name of the output file (default 'line_profile.txt')

        Returns:
        data -- calculated data

        """

        nu, Fline = self.calc_line_profile(nu_min, nu_max, npoints = npoints)

        if include_abs:
            Fabs = self.calc_line_profile(nu_min, nu_max, npoints = npoints, mode = "abs")[-1]

        if include_emit:
            Femit = self.calc_line_profile(nu_min, nu_max, npoints = npoints, mode = "emit")[-1]

        if vs_nu:
            x = nu
            header = "# frequency [Hz]"
        else:
            x = c / nu * 1e8
            header = "# wavelength [Angstrom]"

        data = [nu, Fline]


        header += "   line flux [erg/cm^2/s/Hz]\n"
        if include_abs:
            data.append(Fabs)
            header += "   absorbed flux [erg/cm^2/s/Hz]"
        if include_emit:
            data.append(Femit)
            header += "   emitted flux [erg/cm^2/s/Hz]"

        data = np.array(data)


        f = open(fname, "w")

        f.write("%d spectral points\n" % npoints)
        f.write(header)

        np.savetxt(f, data.T)

        f.close()

        return data

def example():
    """
    Example routine to test the homologous_sphere class
    """

    Lmin = 6.96e10
    Lmax = 6.96e11
    lambda_min = 1200 * 1e-8
    lambda_max = 1230 * 1e-8
    nu_min = c / lambda_max
    nu_max = c / lambda_min
    vmax = 0.01 * c
    vref = 1e8
    ve = 1e8
    tauref = 17.49
    lam0 = 1216.7 * 1e-8

    test = homologous_sphere(rmin = Lmin, rmax = Lmax, vmax = vmax, Ip = 1, tauref = tauref, vref = vref, ve = ve, lam0 = lam0)
    test.show_line_profile(nu_min, nu_max, vs_nu = False)
    test.save_line_profile(nu_min, nu_max, vs_nu = False)

def example_2():
    """
    Example routine to test the homologous_sphere class
    """

    Lmin = 6.96e10
    Lmax = 6.96e11
    lambda_min = 1200 * 1e-8
    lambda_max = 1230 * 1e-8
    nu_min = c / lambda_max
    nu_max = c / lambda_min
    vmax = 0.01 * c
    vref = 1e8
    ve = 1e8
    tauref = 17.49
    lam0 = 1216.7 * 1e-8

    test = homologous_sphere(rmin = Lmin, rmax = Lmax, vmax = vmax, Ip = 1, tauref = tauref, vref = vref, ve = ve, lam0 = lam0)
    test.show_line_profile(nu_min, nu_max, vs_nu = False)
    test.save_line_profile(nu_min, nu_max, vs_nu = False)



if __name__ == "__main__":

    example()
    plt.show()
