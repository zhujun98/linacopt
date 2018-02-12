#!/usr/bin/python

"""
Hold two classes:

- PhaseSpace
    Store and update the phase-space data and the phase-space parameters

- PhaseSpaceParser
    Read the phase-space data from a file


Note: In the calculation (except the canonical emittance), the 
particles are drifted back to the center of the bunch without
considering the collective effects!!!


Author: Jun Zhu

"""
import os

import numpy as np
import pandas as pd

from . import BeamData

V_LIGHT = 299792458
M_E = 9.10938356e-31
Q_E = 1.60217662e-19
CONST_E = M_E*V_LIGHT**2/Q_E 

INF = 1.0e21


class PhaseSpace(BeamData):
    """Store the particle phase-space data and parameters

    Attributes
    ----------
    data: pandas.DataFrame object
        Columns: x (m), y (m), z (m), px (mc), py (mc), pz (mc), t (s, dt)

    particle_file: string
        Name of the particle file.
    particle_type: string
        Type of the particle file.
    _min_pars: int
        Minimum allowed number of particles in the data. The beam
        parameters will not be calculated if the number of particles
        is below this value.

    slice_percent: float
        Percentage of particle for slice parameters calculation.
    cut_halo: None/float
        Percentage of halo to be cut.
    cut_tail: None/float
        Percentage of tail to be cut.
    current_bins: int/'auto'
        No. of bins to calculate the current profile.
    filter_size: int/float
        Standard deviation of the Gaussian kernel of the 1D Gaussian
        filter used for current profile calculation.

    n0: int
        Number of particles in the original file.
    n: int
        Number of particles after applying cut_halo and/or cut_tail.
    charge: int/float
        Charge (C) of the beam after applying cut_halo and/or cut_tail.
    q_norm: int/float/None
        Charge per macro-particle. Only for Impact data.
    p: float
        Average normalized momentum.
    gamma: float
        Average Lorentz factor.
    chirp: float
        Energy chirp (1/m), defined as -<z*delta>/<z^2>.
    Sdelta: float
        Fractional energy spread.
    Sdelta_slice: float
        Fractional slice energy spread.
    St: float
        RMS bunch duration (s).
    dt_slice: float
        Slice full duration (s).
    Sz: float
        RMS bunch length (m).
    I_peak: float
        Peak current (A).
    Sdelta_un: float
        Uncorrelated momentum spread.
    emitx/emity: float
        Normalized horizontal/vertical emittance (m.rad).
    emitx_slice/emity_slice: float
        Normalized horizontal/vertical slice emittance (m.rad).
    Sx/Sy: float
        RMS horizontal/vertical beam size (m).
    betax/betay: float
        Horizontal/vertical beta function (m).
    alphax/alphay: float
        Horizontal/vertical alpha function (m).
    emitx_tr/emity_tr: float
        Normalized horizontal/vertical trace space emittance (m.rad)
    Cx/Cy: float
        Horizontal/vertical displacement (m).
    Cxp/Cyp: float
        Horizontal/vertical divergence (rad).
    Ct: float
        Average timing (s).
    """
    def __init__(self, particle_file, particle_type, charge=None, q_norm=None,
                 slice_percent=0.1, cut_halo=None, cut_tail=None, rotate=None,
                 current_bins='auto', filter_size=1, min_pars=5, opt=False,
                 slice_with_peak_current=True):
        """
        Parameters
        ----------
        particle_file: string
            Name of the particle file.
        particle_type: string
            Type of the particle file.
        charge: int/float
            Charge of the beam (in C). Only for Impact data.
            Ignored if q_norm is given.
        q_norm: int/float/None
            Charge per macro-particle. Only for Impact data.
        slice_percent: float
            Percent of the slice bunch length to the total bunch length.
        cut_halo: None/float
            Percentage of particles to be removed based on their
            transverse distance to the bunch centroid. Applied
            before tail cutting.
        cut_tail: None/float
            Percentage of particles to be removed in the tail.
        rotate: None/float
            Apply rotation (in degree) to the phasespace.
        current_bins: int/'auto'
            No. of bins to calculate the current profile.
        filter_size: int/float
            Standard deviation of the Gaussian kernel of the 1D Gaussian
            filter used for current profile calculation.
        min_pars: int
            Minimum allowed number of particles in the data.
        opt: Boolean
            True for the initialization of the fit-points in linac_opt.
            Since there is no output, an error will occur if the update
            method is called. Default is False.
        use_slice_with_peak_current: Boolean
            True for calculating slice properties of the slice with peak
            current; False for calculating slice properties of the slice
            in the center of the bunch.
        """
        self.particle_file = os.path.join(os.getcwd(), particle_file)
        self._min_pars = min_pars

        super().__init__(particle_type)

        self.charge = charge
        self.q_norm = q_norm
        self.current_dist = None

        if type(slice_percent) in (int, float) and 0.0 < slice_percent < 1.0:
            self.slice_percent = slice_percent
        else:
            self.slice_percent = 0.1
        self.current_bins = current_bins
        self.filter_size = filter_size

        self.cut_halo = cut_halo
        self.cut_tail = cut_tail
        self.rotate = rotate

        self.data = None

        self.n0 = None
        self.n = None

        self.p = None
        self.gamma = None
        self.chirp = None
        self.Sdelta = None
        self.Sdelta_slice = None
        self.St = None
        self.dt_slice = None
        self.Sz = None
        self.I_peak = None
        self.Sdelta_un = None

        self.emitx = None
        self.emitx_slice = None
        self.Sx = None
        self.betax = None
        self.alphax = None
        self.emitx_tr = None

        self.emity = None
        self.emity_slice = None
        self.Sy = None
        self.betay = None
        self.alphay = None
        self.emity_tr = None

        self.Cx = None
        self.Cxp = None
        self.Cy = None
        self.Cyp = None
        self.Ct = None

        self.slice_with_peak_current = slice_with_peak_current

        if not opt:
            self.update()

    def update(self):
        """Read the phase-space and calculate the beam parameters"""
        self.data = None
        parser = PhaseSpaceParser()
        if self.particle_type == 'astra':
            self.data, self.charge = parser.astra_parser(self.particle_file)
        elif self.particle_type == 'impact':
            self.data = parser.impact_parser(self.particle_file)

        if isinstance(self.rotate, (float, int)):
            self._rotate()

        self.n0 = len(self.data)
        self.n = self.n0

        if self.particle_type == 'impact':
            if type(self.q_norm) in [float, int]:
                self.charge = self.q_norm*self.n0
            if self.charge is None:
                self.charge = 0.0

        # Cut the halo of the bunch
        if isinstance(self.cut_halo, float) and 0.0 < self.cut_halo < 1.0:
            self.n = int(self.n*(1 - self.cut_halo))
            self.charge *= float(self.n) / self.n0

            self.data['r'] = np.sqrt(self.data['x'] ** 2 + self.data['y'] ** 2)
            self.data = self.data.reindex(
                self.data['r'].sort_values(ascending=True).index)
            self.data = self.data[:self.n]
            del self.data['r']

        # Cut the tail of the bunch
        if isinstance(self.cut_tail, float) and 0.0 < self.cut_tail < 1.0:
            self.n = int(self.n*(1 - self.cut_tail))
            self.charge *= float(self.n)/self.n0
            # User median() here to deal with extreme outliers
            self.data['t'] -= self.data['t'].median()
            self.data = self.data.reindex(
                self.data['t'].abs().sort_values(ascending=True).index)
            self.data = self.data[:self.n]
            self.data['t'] -= self.data['t'].mean()

        # Too few particles may cause error during the following
        # calculation, e.g. negative value in sqrt.
        if self.n < self._min_pars:
            raise ValueError("Too few particles ({}) in the data!".
                             format(self.n))

        p = np.sqrt(self.data['pz']**2 + self.data['px']**2 + self.data['py']**2)

        p_ave = p.mean()
        dp = (p - p_ave) / p_ave
        dz = self.data['z'] - self.data['z'].mean()

        self.p = p_ave
        self.gamma = np.sqrt(p_ave ** 2 + 1)
        self.chirp = -1 * dp.cov(dz) / dz.var(ddof=0)
        self.Sdelta = p.std(ddof=0) / p_ave
        self.St = self.data['t'].std(ddof=0)
        self.Sz = self.data['z'].std(ddof=0)

        # The current profile calculation is included here but not in
        # the plot function. If it is included in the plot function, the
        # printout parameters might differ from the plot, which could
        # cause confusion.
        counts, edges = np.histogram(self.data['t'], bins=self.current_bins)
        step_size = edges[1] - edges[0]
        centers = edges[:-1] + step_size / 2
        current = counts / float(len(self.data)) * self.charge / step_size
        current = self.gaussian_filter1d(current, sigma=self.filter_size)
        self.I_peak = current.max()
        self.current_dist = [centers, current]

        self.emitx = self._canonical_emit(self.data['x'], self.data['px'])

        self.Sx, self.betax, self.alphax, self.emitx_tr \
            = self._twiss(self.data['x'], dz, self.data['px'],
                          self.data['pz'], self.gamma)

        self.emity = self._canonical_emit(self.data['y'], self.data['py'])

        self.Sy, self.betay, self.alphay, self.emity_tr \
            = self._twiss(self.data['y'], dz, self.data['py'],
                          self.data['pz'], self.gamma)

        self.Cx = self.data['x'].mean()
        self.Cy = self.data['y'].mean()
        self.Cxp = (self.data['px'] / self.data['pz']).mean()
        self.Cyp = (self.data['py'] / self.data['pz']).mean()
        self.Ct = self.data['t'].mean()

        # Calculate the slice parameters
        sorted_data = self.data.reindex(
            self.data['t'].abs().sort_values(ascending=True).index)

        if self.slice_with_peak_current is True:
            Ct_slice = centers[np.argmax(current)]
        else:
            Ct_slice = self.Ct

        dt_slice = 4*self.St*self.slice_percent  # assume 4-sigma full bunch length
        slice_data = sorted_data[(sorted_data.t > Ct_slice - dt_slice/2) &
                                 (sorted_data.t < Ct_slice + dt_slice/2)]

        if len(slice_data) < self._min_pars:
            raise ValueError("Too few particles ({}) in the slice data!".
                             format(len(slice_data)))

        p_slice = np.sqrt(slice_data['pz']**2 + slice_data['px']**2
                          + slice_data['py']**2)

        self.emitx_slice = self._canonical_emit(slice_data.x, slice_data.px)
        self.emity_slice = self._canonical_emit(slice_data.y, slice_data.py)
        self.Sdelta_slice = p_slice.std(ddof=0) / p_slice.mean()
        self.dt_slice = slice_data.t.max() - slice_data.t.min()

        # The output will be different from the output by msddsplot
        # because the slightly different No of particles sliced. It
        # affects the correlation calculation since the value is
        # already very close to 1.
        self.Sdelta_un = self.Sdelta_slice * np.sqrt(
            1 - (slice_data['t'].corr(p_slice)) ** 2)

    def _rotate(self):
        """Rotate the phasespace

        Provided by F. Mayer
        """
        theta = self.rotate * np.pi / 180.0  # Convert to rad

        # x(m), px(mc), y(m), py(mc), t(s), p(mc), z(m).

        cos_theta = np.cos(theta)
        sin_theta = np.sin(theta)

        def transformation(r, cm):
            x = r[0]
            y = r[1]
            z = r[2]
            cx = cm[0]
            cy = cm[1]
            cz = cm[2]

            x_new = cx - cx*cos_theta + x*cos_theta - cz*sin_theta + z*sin_theta
            y_new = y
            z_new = cz - cz*cos_theta + z*cos_theta + cx*sin_theta - x*sin_theta

            return [x_new, y_new, z_new]

        pos = [self.data['x'], self.data['y'], self.data['z']]
        cm_pos = [np.mean(self.data['x']), np.mean(self.data['y']),
                  np.mean(self.data['z'])]
        mom = [self.data['px'], self.data['py'], self.data['pz']]
        cm_mom = [np.mean(self.data['px']), np.mean(self.data['py']),
                  np.mean(self.data['pz'])]

        [self.data['x'], self.data['y'], self.data['z']] = \
            transformation(pos, cm_pos)
        [self.data['px'], self.data['py'], self.data['pz']] = \
            transformation(mom, cm_mom)

    @staticmethod
    def _canonical_emit(x, px):
        """ Calculate the canonical emittance.

        Parameters
        ----------
        x: pandas.Series object
            Position coordinates
        px: pandas.Series object
            Momentum coordinates

        Return
        ------
        Normalized canonical emittance.
        """
        x_ave = x.mean()
        px_ave = px.mean()
        x2 = x.var(ddof=0)
        px2 = px.var(ddof=0)
        xpx = ((x - x_ave) * (px - px_ave)).mean()

        return np.sqrt(x2 * px2 - xpx ** 2)

    @staticmethod
    def _twiss(x, dz, px, pz, gamma):
        """ Calculate the Twiss parameters

        The particles are drifted back to the centroid of the bunch!

        Parameters
        ----------
        x: pandas.Series object
            Position coordinates
        dz: pandas.Series object
            Longitudinal distance to the bunch centre.
        px: pandas.Series object
            Momentum coordinates
        pz: pandas.Series object
            Longitudinal momentum
        gamma: float
            Average Lorentz factor of the bunch

        Returns
        -------
        sigma_x: float
            RMS transverse beam size.
        betax: float
            Beta function.
        alphax: float
            Alpha function.
        emitnx: float
            Normalized trace-space emittance
        """
        beta = np.sqrt(1 - 1 / gamma ** 2)

        xp = px / pz
        x_new = x - dz * xp

        x_ave = x_new.mean()
        xp_ave = xp.mean()
        x2 = x_new.var(ddof=0)
        xp2 = xp.var(ddof=0)
        xxp = ((x_new - x_ave) * (xp - xp_ave)).mean()

        emitx = np.sqrt(x2 * xp2 - xxp ** 2)
        emitnx = emitx * beta * gamma
        sigma_x = np.sqrt(x2)
        betax = x2 / emitx
        alphax = -1 * xxp / emitx

        return [sigma_x, betax, alphax, emitnx]

    @staticmethod
    def gaussian_filter1d(x, sigma):
        """One-dimensional Gaussian filter.

        Parameters
        ----------
        x: array_like
            Input array for filter
        sigma: int/float
            Standard deviation for Gaussian kernel.

        Returns
        -------
        Filtered x.
        """
        sd = float(sigma)
        # make the radius of the filter equal to truncate standard deviations
        lw = int(4.0*sd + 0.5)
        weights = [0.0]*(2*lw + 1)
        weights[lw] = 1.0
        sum_ = 1.0
        sd *= sd
        # calculate the kernel:
        for ii in range(1, lw + 1):
            tmp = np.exp(-0.5 * float(ii * ii) / sd)
            weights[lw + ii] = tmp
            weights[lw - ii] = tmp
            sum_ += 2.0 * tmp

        weights /= sum_

        return np.convolve(x, weights, 'same')

    def output_params(self, file_name='params.out'):
        """Print the beam parameters into a file."""
        output = os.path.join(os.path.dirname(self.particle_file), file_name)
        with open(output, 'w') as f:
            print(self, file=f)
            
        print(("Saved parameters at {}".format(file_name)))

    def __str__(self):
        """"""
        text = "cutTail = {}, cutHalo = {}, rotate = {}\n\n". \
            format(self.cut_tail, self.cut_halo, self.rotate)
        text += "{:16}    {:16}    {:16}    {:16}\n". \
            format('n', 'charge (C)', 'p', 'I_peak (A)')
        text += "{:16.4e}    {:16.4e}    {:16.4e}    {:16.4e}\n\n". \
            format(self.n, self.charge, self.p, self.I_peak)
        text += "{:16}    {:16}    {:16}    {:16}\n". \
            format('emitx (m)', 'emity (m)', 'Sx (m)', 'Sy (m)')
        text += "{:16.4e}    {:16.4e}    {:16.4e}    {:16.4e}\n\n". \
            format(self.emitx, self.emity, self.Sx, self.Sy)
        text += "{:16}    {:16}    {:16}    {:16}\n". \
            format('betax (m)', 'betay (m)', 'alphax', 'alphay')
        text += "{:16.4e}    {:16.4e}    {:16.4e}    {:16.4e}\n\n". \
            format(self.betax, self.betay, self.alphax, self.alphay)
        text += "{:16}    {:16}    {:16}    {:16}\n". \
            format('St (s)', 'Sdelta', 'chirp (1/m)', 'Ct (s)')
        text += "{:16.4e}    {:16.4e}    {:16.4e}    {:16.4e}\n\n". \
            format(self.St, self.Sdelta, self.chirp, self.Ct)
        text += "{:16}    {:16}    {:16}    {:16}\n". \
            format('emitx_slice (m)', 'emity_slice (m)', 'Sdelta_slice', 'dt_slice (s)')
        text += "{:16.4e}    {:16.4e}    {:16.4e}    {:16.4e}\n\n". \
            format(self.emitx_slice, self.emity_slice, self.Sdelta_slice, self.dt_slice)
        text += "{:16}    {:16}    {:16}    {:16}\n". \
            format('Cx (m)', 'Cy (m)', 'Cxp (rad)', 'Cyp (rad)')
        text += "{:16.4e}    {:16.4e}    {:16.4e}    {:16.4e}\n\n". \
            format(self.Cx, self.Cy, self.Cxp, self.Cyp)
        text += "{:16}    {:16}    {:16}\n". \
            format('emitx_tr (m)', 'emity_tr (m)', 'Sdelta_un')
        text += "{:16.4e}    {:16.4e}    {:16.4e}\n". \
            format(self.emitx_tr, self.emity_tr, self.Sdelta_un)

        return text


class PhaseSpaceParser(object):
    """Read the phase-space data from file

    The returned phase-space data contains the following columns:
        x (m), px (mc), y (m), py (mc), z(m), pz (mc), t (s).
    """
    @staticmethod
    def astra_parser(particle_file):
        """Parse the ASTRA particle file.

        Parameter
        ---------
        particle_file: string
            Name of the particle file

        Return
        ------
        data: pandas.DataFrame
            Phase-space data.
        charge: float
            Charge (C) of the bunch.
        """
        # Units: m, m, m, eV/c, eV/c, eV/c, ns, nC, NA, NA
        col_names = ['x', 'y', 'z', 'px', 'py', 'pz',
                     't', 'charge', 'index', 'flag']

        data = pd.read_csv(particle_file, delim_whitespace=True,
                           names=col_names)

        pz_ref = data['pz'].iloc[0]
        data.ix[0, 'pz'] = 0.0
        data['pz'] += pz_ref

        data['px'] /= M_E * V_LIGHT ** 2 / Q_E
        data['py'] /= M_E * V_LIGHT ** 2 / Q_E
        data['pz'] /= M_E * V_LIGHT ** 2 / Q_E

        # ix will first try to act like loc to find the index label.
        # If the index label is not found, it will add an index label
        # (new row). Therefore, the reference particle must be used
        # before removing lost particles since the reference particle
        # could be removed.
        z_ref = data['z'].iloc[0]
        data.ix[0, 'z'] = 0.0
        data['z'] += z_ref

        # remove lost particles
        data = data[data['flag'].isin([3, 5])]

        data['p'] = np.sqrt(data['px'] ** 2 + data['py'] ** 2 + data['pz'] ** 2)

        # At this step, the timing can be used for timing jitter study.
        data['t'] = data['t'].iloc[0]/1e9 - (data['z'] - z_ref)\
            /(V_LIGHT * data['pz']/np.sqrt(data['p']**2 + 1))

        # The bunch is centered for the convenience of the longitudinal
        # phase-space plot.
        data['t'] = data['t'] - data['t'].mean()

        charge = -1e-9 * data['charge'].sum()

        data.drop(['charge', 'index', 'flag'], inplace=True, axis=1)

        return data, charge

    @staticmethod
    def impact_parser(particle_file):
        """Parse the IMPACT-T/Z particle file.

        Parameter
        ---------
        particle_file: string
            Name of the particle file

        Return
        ------
        data: pandas.DataFrame
            Phase-space data.
        """
        # Units: m, /mc, m, /mc, m, /mc
        col_names = ['x', 'px', 'y', 'py', 'z', 'pz']

        data = pd.read_csv(particle_file, delim_whitespace=True, names=col_names)

        # Drop the first row if the input file is 'partcl.data'.
        data.dropna(inplace=True)

        data['p'] = \
            np.sqrt(data['px'] ** 2 + data['py'] ** 2 + data['pz'] ** 2)

        data['t'] = -(data['z'] - data['z'].mean()) \
                    / (V_LIGHT * data['pz'] / np.sqrt(data['p'] ** 2 + 1))

        # data.drop(['p'], inplace=True, axis=1)

        return data


if __name__ == "__main__":
    # ps_astra = PhaseSpace('examples/plots/injector.0400.001', 'astra')
    ps_astra = PhaseSpace('examples/astra_advanced/injector.0620.001', 'astra')
    print('-'*80 + "\nParameters for {}:\n".format(ps_astra.particle_file))
    print(ps_astra)
    ps_astra.output_params()

    # ps_impact = PhaseSpace('examples/plots/fort.140', 'impact',
    #                        charge=17.7e-12, q_norm=None, cut_tail=0.0, rotate=0.0,
    #                        cut_halo=0.0, current_bins='auto', filter_size=1,
    #                        slice_percent=0.1, min_pars=10,
    #                        slice_with_peak_current=True)
    # ps_impact.update()
    # print('-'*80 + "\nParameters for {}:\n".format(ps_impact.particle_file))
    # print(ps_impact)
