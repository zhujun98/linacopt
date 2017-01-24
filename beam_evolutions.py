#!/usr/bin/python

"""
beam_evolutions

Hold three classes:

    - BeamEvolution:
        Store and update the beam evolution and its statistics

    - BeamEvolutionParser
        Read the beam evolution data from files.

    - Stats
        Store and update the statistics of a single variable

Important note:
_______________
The calculated Twiss parameters are based on the canonical coordinates
of the beam. Therefore, when the beam has a large energy spread or a
big divergence angle, the result will be quite different from the
true Twiss parameters. However, one can set TR_EmitS = .T in ASTRA to
get the very-close Twiss parameters even in extreme conditions.

Tested on:
__________
Ubuntu 14.04
Ubuntu 16.04

Developer:
__________
Jun Zhu

Last modified on 22/01/2017

"""
import os

import numpy as np
import pandas as pd

from linacopt_data import LinacOptData


V_LIGHT = 299792458
M_E = 9.10938356e-31
Q_E = 1.60217662e-19
CONST_E = M_E*V_LIGHT**2/Q_E

INF = 1.0e21


class BeamEvolution(LinacOptData):
    """Store the beam evolution and its statistics

    Attributes
    ----------
    data: pandas.DataFrame object
        z (m), gamma, SdE (eV), Sx (m), Sy (m), Sz (m),
        emitx (m.rad), emity (m.rad), emitz (m.rad),
        emitx_tr (m.rad), emity_tr (m.rad),
        betax (m), betay (m), alphax, alphay

    z: Stats object
        Longitudinal position (m).
    gamma: Stats object
        Lorentz factor.
    SdE: Stats object
        RMS energy spread (eV).
    Sx/Sy/Sz: Stats objects
        RMS bunch sizes (m).
    betax/betay: Stats objects
        Beta functions (m).
    alphax/alphay: Stats objects
        Alpha functions.
    emitx/emity: Stats objects
        Normalized canonical emittance (m.rad).
    emitx_tr/emity_tr: Stats objects
        Normalized trace-space emittance (m.rad).
    """
    def __init__(self, root_name, particle_type, z_lim=None, opt=False):
        """Initialize BeamStats object

        Parameters
        ----------
        root_name: string
            The root name of the output files. For Impact-T files,
            root_name will be set to 'fort' if not given.
        particle_type: string
            Type of the particle file.
        z_lim: scalar/tuple
            If None, passed as (-INF, INF)
            if scalar, being passed as (left, INF)
            if tuple, the first two elements being passed as (left, right)
        opt: Boolean
            True for the initialization of the fit-points in linac_opt.
            Since there is no output, an error will occur if the update
            method is called. Default is False.
        """
        self.root_name = os.path.join(os.getcwd(), root_name)

        super(BeamEvolution, self).__init__(particle_type)

        if z_lim is None:
            self.z_lim = (-1*INF, INF)
        else:
            try:
                self.z_lim = tuple(z_lim)[0:2]
            except TypeError:
                self.z_lim = (z_lim, INF)

        self.data = None

        self.z = Stats()
        self.gamma = Stats()
        self.SdE = Stats()
        self.Sx = Stats()
        self.Sy = Stats()
        self.Sz = Stats()
        self.betax = Stats()
        self.betay = Stats()
        self.alphax = Stats()
        self.alphay = Stats()
        self.emitx = Stats()
        self.emity = Stats()
        self.emitx_tr = Stats()
        self.emity_tr = Stats()

        if not opt:
            self.update()

    def update(self):
        """Update attributes of attributes"""
        parser = BeamEvolutionParser()
        if self.particle_type == 'astra':
            self.data = parser.astra_parser(self.root_name)
        elif self.particle_type == 'impact':
            self.data = parser.impact_parser(self.root_name)
        else:
            raise ValueError(
                "The particle_type must be either 'astra' or 'impact'!\n")

        self.slice_data()

        for key, value in self.__dict__.iteritems():
            if isinstance(value, Stats):
                value.update(self.data[key])

    def slice_data(self):
        """Slice data in the range of self.z_lim"""
        if self.z_lim[0] >= self.data['z'].max() or \
                self.z_lim[1] <= self.data['z'].min():
            raise ValueError("\nOut of range: z_lim!")

        i_min = 0
        i_max = len(self.data['z'])
        for i in range(i_max):
            if self.z_lim[0] <= self.data['z'][i]:
                i_min = i
                break

        for i in range(i_max):
            if self.z_lim[1] < self.data['z'][i]:
                i_max = i - 1
                break

        self.data = self.data.ix[i_min:i_max]

    def __str__(self):
        """Print output"""
        text = '- z (m)\n'
        text += str(self.__getattribute__('z'))
        text += '- gamma\n'
        text += str(self.__getattribute__('gamma'))
        text += '- SdE (eV)\n'
        text += str(self.__getattribute__('SdE'))
        text += '- Sx (m)\n'
        text += str(self.__getattribute__('Sx'))
        text += '- Sy (m)\n'
        text += str(self.__getattribute__('Sy'))
        text += '- Sz (m)\n'
        text += str(self.__getattribute__('Sz'))
        text += '- betax (m)\n'
        text += str(self.__getattribute__('betax'))
        text += '- betay (m)\n'
        text += str(self.__getattribute__('betay'))
        text += '- alphax (m)\n'
        text += str(self.__getattribute__('alphax'))
        text += '- alphay (m)\n'
        text += str(self.__getattribute__('alphay'))
        text += '- emitx (m.rad)\n'
        text += str(self.__getattribute__('emitx'))
        text += '- emity (m.rad)\n'
        text += str(self.__getattribute__('emity'))
        text += '- emitx_tr (m.rad)\n'
        text += str(self.__getattribute__('emitx_tr'))
        text += '- emity_tr (m.rad)\n'
        text += str(self.__getattribute__('emity_tr'))

        return text


class BeamEvolutionParser(object):
    """Read the line data from files

    The lina data is a pandas.DataFrame with the following columns:
        z (m), gamma, SdE (eV), Sx (m), Sy (m), Sz (m),
        emitx (m.rad), emity (m.rad), emitz (m.rad),
        emitx_tr (m.rad), emity_tr (m.rad),
        betax (m), betay (m), alphax, alphay
    """
    @staticmethod
    def impact_parser(root_name):
        """Parse the IMPACT-T/Z particle file.

        Columns in the input files
        --------------------------
        xdata:
            t (s), z (m), Cx (m), Sx (m), px (/mc), Spx (/mc),
            twiss (m), emitx (m).
        ydata:
            t (s), z (m), Cy (m), Sy (m), py (/mc), Spy (/mc),
            twiss (m), emity (m).
        zdata:.
            t (s), z (m), Sz (m), pz (/mc), Spz (/mc),
            twiss (m), emitz (m).

        Note: twiss = -<x - <x>><px - <px>>
        """
        if root_name is None:
            root_name == 'fort'

        x_file = root_name + '.24'
        y_file = root_name + '.25'
        z_file = root_name + '.26'

        xdata = pd.read_csv(
            x_file, delim_whitespace=True,
            names=['t', 'z', 'Cx', 'Sx', 'px', 'Spx', 'x_px', 'emitx'])
        ydata = pd.read_csv(
            y_file, delim_whitespace=True,
            names=['t', 'z', 'Cy', 'Sy', 'py', 'Spy', 'y_py', 'emity'])
        zdata = pd.read_csv(
            z_file, delim_whitespace=True,
            names=['t', 'z', 'Sz', 'pz', 'Spz', 'z_pz', 'emitz'])

        line_data = pd.DataFrame()

        line_data['z'] = xdata['z']
        # line_data['t'] = xdata['t']

        line_data['gamma'] = np.sqrt(xdata['px'].pow(2) + ydata['py'].pow(2) +
                                     zdata['pz'].pow(2) + 1)
        line_data['SdE'] = np.sqrt(xdata['Spx'].pow(2) + ydata['Spy'].pow(2) +
                                   zdata['Spz'].pow(2))*CONST_E

        line_data['Sx'] = xdata['Sx']
        line_data['Sy'] = ydata['Sy']
        line_data['Sz'] = zdata['Sz']

        line_data['emitx'] = xdata['emitx']
        line_data['emity'] = ydata['emity']
        line_data['emitz'] = zdata['emitz']

        line_data['emitx_tr'] = xdata['emitx']
        line_data['emity_tr'] = ydata['emity']
        # line_data['emitz_tr'] = zdata['emitz']

        gamma_beta = line_data['gamma']*np.sqrt(1 - 1/line_data['gamma'].pow(2))

        line_data['betax'] = line_data['Sx'].pow(2)*gamma_beta/line_data['emitx_tr']
        line_data['betay'] = line_data['Sy'].pow(2)*gamma_beta/line_data['emity_tr']

        line_data['alphax'] = xdata['x_px']/line_data['emitx_tr']
        line_data['alphay'] = ydata['y_py']/line_data['emity_tr']

        return line_data

    @staticmethod
    def astra_parser(root_name):
        """Parse the ASTRA particle file.

        Columns in the input files
        --------------------------
        xdata:
            z (m), t (ns), Cx (mm), Sx (mm), Sxp (mm), emitx (um),
            x_xp (um).
        ydata:
            z (m), t (ns), Cy (mm), Sy (mm), Syp (mm), emity (um),
            y_yp (um).
        zdata:
            z (m), t (ns), Ek (MeV), Sz (mm), SdE (keV), emitz (um),
            z_dE (um).

        Note: x_xp = <x.*xp>/sqrt(<x^2>)/<pz>,
              y_yp = <y.*yp>/sqrt(<y^2>)/<pz>,
              z_de = <z.*dE>/sqrt(<z^2>),
              Sxp = Spx/<pz>, this is not Sxp, so I cannot calculate
              the trace-space emittance!!!
              Syp = Spy/<pz>.

              Even thought the trace-space emittance is read from the
              TRemit file, there is still a little difference between
              the correct value since we do not know the real <x.*xp>
        """
        if root_name is None:
            raise ValueError("\nroot_name of the output files is not given!")

        x_file = root_name + '.Xemit.001'
        y_file = root_name + '.Yemit.001'
        z_file = root_name + '.Zemit.001'
        emit_tr_file = root_name + '.TRemit.001'

        line_data = pd.DataFrame()

        xdata = pd.read_csv(
            x_file, delim_whitespace=True,
            names=['z', 't', 'Cx', 'Sx', 'Sxp', 'emitx', 'x_xp'])
        ydata = pd.read_csv(
            y_file, delim_whitespace=True,
            names=['z', 't', 'Cy', 'Sy', 'Syp', 'emity', 'y_yp'])
        zdata = pd.read_csv(
            z_file, delim_whitespace=True,
            names=['z', 't', 'Ek', 'Sz', 'SdE', 'emitz', 'z_dE'])

        # ASTRA will not output .TRemit file by default
        try:
            emit_tr_data = pd.read_csv(
                emit_tr_file, delim_whitespace=True,
                names=['z', 't', 'emitx_tr', 'emity_tr', 'emitz_tr'])
            line_data['emitx_tr'] = emit_tr_data['emitx_tr']*1.0e-6
            line_data['emity_tr'] = emit_tr_data['emity_tr']*1.0e-6
            # line_data['emitz_tr'] = emit_tr_data['emitz_tr']*1.0e-6
        except IOError:
            line_data['emitx_tr'] = xdata['emitx']*1.0e-6
            line_data['emity_tr'] = ydata['emity']*1.0e-6
            # line_data['emitz_tr'] = zdata['emitz']*1.0e-6

        line_data['z'] = xdata['z']
        # line_data['t'] = xdata['t']*1.0e-9

        line_data['gamma'] = zdata['Ek']*1.0e6/CONST_E + 1
        line_data['SdE'] = zdata['SdE']*1.0e3

        line_data['Sx'] = xdata['Sx']*1.0e-3
        line_data['Sy'] = ydata['Sy']*1.0e-3
        line_data['Sz'] = zdata['Sz']*1.0e-3

        # line_data['Sxp'] = xdata['Sxp']*1.0e-3
        # line_data['Syp'] = ydata['Syp']*1.0e-3

        line_data['emitx'] = xdata['emitx']*1.0e-6
        line_data['emity'] = ydata['emity']*1.0e-6
        line_data['emitz'] = zdata['emitz']*1.0e-6

        gamma_beta = line_data['gamma']*np.sqrt(1 - 1/line_data['gamma'].pow(2))

        line_data['betax'] = line_data['Sx'].pow(2)*gamma_beta/line_data['emitx_tr']
        line_data['betay'] = line_data['Sy'].pow(2)*gamma_beta/line_data['emity_tr']

        x_xp = xdata['Sx']*xdata['x_xp']*1.0e-6
        y_yp = ydata['Sy']*ydata['y_yp']*1.0e-6
        line_data['alphax'] = -x_xp*gamma_beta/line_data['emitx_tr']
        line_data['alphay'] = -y_yp*gamma_beta/line_data['emity_tr']

        return line_data


class Stats(object):
    """Store the statistic values of an array-like object.

    Attributes
    ----------
    start: float
        First value.
    end: float
        Last value.
    max: float
        Maximum value.
    min: float
        Minimum value.
    ave: float
        Average value.
    std: float
        Standard deviation.
    """
    def __init__(self):
        """"""
        self.start = None
        self.end = None
        self.max = None
        self.min = None
        self.ave = None
        self.std = None

    def __str__(self):
        """"""
        text = "{:12}    {:12}    {:12}    {:12}    {:12}    {:12}\n".\
            format('start', 'end', 'minimum', 'maximum', 'average', 'std')

        text += "{:12.4e}    {:12.4e}    {:12.4e}    {:12.4e}    {:12.4e}    " \
                "{:12.4e}\n\n".format(self.start, self.end, self.min, self.max,
                                      self.ave, self.std)

        return text

    def update(self, data):
        """Update attributes

        Parameter
        ---------
        data: array-like
            Input data.
        """
        data = np.asarray(data)
        if data.ndim > 1:
            raise ValueError("One-dimensional array is foreseen!")

        self.start = data[0]
        self.end = data[-1]
        self.max = data.max()
        self.min = data.min()
        self.ave = data.mean()
        self.std = data.std(ddof=0)


if __name__ == "__main__":
    # Test
    ld_astra = BeamEvolution('astra_test/injector', 'astra', z_lim=(4, 6))
    print '-'*80 + "\nParameters for test ASTRA files"
    print ld_astra

    ld_impact = BeamEvolution('impact_test/fort', 'impact')
    print '-'*80 + "\nParameters for test IMPACT files"
    print ld_impact
