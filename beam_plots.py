#!/usr/bin/python
"""
Hold two classes:

- PhaseSpacePlot:
    Plot the beam phase-space.

- LinePlot:
    Plot beam parameter evolutions along the beamline.


Author: Jun Zhu


TODO:
- Include lattice in the line plot.
- matplotlib.rcParams in different systems.
"""
import os
import re
import random

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.ndimage.filters import gaussian_filter

from beam_parameters import PhaseSpace
from beam_evolutions import BeamEvolution


LABEL_FONT_SIZE = 26
TITLE_FONT_SIZE = 18
TICK_FONT_SIZE = 20
LEGEND_FONT_SIZE = 18
LABEL_PAD = 8
TICK_PAD = 8
MAX_LOCATOR = 7
AX_MARGIN = 0.05


class PhaseSpacePlot(PhaseSpace):
    """Plot the beam phase-space

    The purpose of creating another class is that people do not have to
    install matplotlib for the optimizer.
    """
    def plot(self, x, y, special=True, x_unit=None, y_unit=None, y1_unit=None,
             x_lim=None, y_lim=None, density_plot=True, marker_size=2, 
             marker_color='dodgerblue', alpha=1.0, bins_2d=500, sigma_2d=5, 
             sample=20000, figsize=None, output='', dpi=300):
        """Plot different phase-spaces

        Parameters
        ----------
        x: string
            Variable in x-axis.
        y: string
            Variable in y-axis.
        special: Boolean
            True for turning on special effects of some plots.
        x_unit: string
            Units of y1 axis.
        y_unit: string
            Units of y1 axis.
        y1_unit: string
            Units of y1 axis.
        x_lim -> tuple
            Limits for the x axis as (left, right)
        y_lim -> tuple
            Limits for the y axis as (left, right)
        density_plot: Boolean
            True for colorful density plot.
        marker_size: int
            Size of the point in the scatter plot.
        marker_color: string
            Color of markers for non-density plot.
        alpha: float, [0, 1]
            Alpha value (transparency).
        bins_2d: int or [int, int]
            No. of bins used in numpy.histogram2d.
        sigma_2d: int/float
            Standard deviation of Gaussian kernel of the Gaussian filter.
        sample: non-negative int/float
            If sample < 1.0, sample by fraction else sample by count
            (round to integer).
        figsize: None/tuple
            Size of the figure (width, height).
        output: string
            Name of the output file.
        dpi: int
            dpi of the png output.
        """
        options = ['x', 'y', 'z', 'xp', 'yp', 't', 'p']
        if x not in options or y not in options:
            raise ValueError("Valid options are: {}".format(options))

        phasespace = self.data
        if phasespace is None:
            raise ValueError("NoneType phase-space data!")

        if figsize is None:
            figsize = (8, 6)
            if density_plot is not True:
                figsize = (6.8, 6)
        fig, ax = plt.subplots(figsize=figsize)

        ax.margins(AX_MARGIN)

        if x_unit is None:
            x_unit = default_unit(x)
        if y_unit is None:
            y_unit = default_unit(y)

        x_unit_label, x_scale = unit_scale(x_unit)
        y_unit_label, y_scale = unit_scale(y_unit)

        x_sample, y_sample, density_color, i_sample = \
            self.sample_data(phasespace[self.get_column(x)], phasespace[self.get_column(y)],
                             bins=bins_2d, sigma=sigma_2d, sample=sample)

        if x in ['xp', 'yp']:
            x_sample /= phasespace.pz.iloc[i_sample]
        if y in ['xp', 'yp']:
            y_sample /= phasespace.pz.iloc[i_sample]

        x_symmetric = False
        y_symmetric = False
        if x in ('x', 'xp'):
            x_symmetric = True
        if y in ('y', 'yp'):
            y_symmetric = True
        ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=MAX_LOCATOR,
                                                      symmetric=x_symmetric))
        ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=MAX_LOCATOR,
                                                      symmetric=y_symmetric))

        if density_plot is False:
            ax.scatter(x_sample*x_scale, y_sample*y_scale, alpha=alpha,
                       c=marker_color, edgecolor='', s=marker_size)
        else:
            cb = ax.scatter(x_sample*x_scale, y_sample*y_scale, c=density_color,
                            alpha=alpha, cmap='jet', edgecolor='', s=marker_size)

            if (x, y) == ('t', 'p') and special is True:
                cbaxes = fig.add_axes([0.75, 0.07, 0.2, 0.02])
                cbar = plt.colorbar(cb, orientation='horizontal', cax=cbaxes)
            else:
                cbar = plt.colorbar(cb, shrink=0.5)
            cbar.set_ticks(np.arange(0, 1.01, 0.2))
            cbar.ax.tick_params(labelsize=14)

        ax.set_xlabel(name2label(x) + ' ' + x_unit_label,
                      fontsize=LABEL_FONT_SIZE, labelpad=LABEL_PAD)
        ax.set_ylabel(name2label(y) + ' ' + y_unit_label,
                      fontsize=LABEL_FONT_SIZE, labelpad=LABEL_PAD)
        ax.tick_params(labelsize=TICK_FONT_SIZE, pad=TICK_PAD)

        if (x, y) == ('x', 'xp') and special is True:
            plt.title(r'$\varepsilon_x$ = %s $\mu$m' %
                      float("%.2g" % (self.emitx*1e6)),
                      fontsize=TITLE_FONT_SIZE, y=1.02)

        if (x, y) == ('y', 'yp') and special is True:
            plt.title(r'$\varepsilon_y$ = %s $\mu$m' %
                      float("%.2g" % (self.emity*1e6)),
                      fontsize=TITLE_FONT_SIZE, y=1.02)

        if (x, y) == ('t', 'p') and special is True:
            if not y1_unit:
                y1_unit = 'A'

            ax1 = ax.twinx()
            ax1.margins(AX_MARGIN)
            ax1.yaxis.set_major_locator(ticker.MaxNLocator(nbins=MAX_LOCATOR))

            y1_unit_label, y1_scale = unit_scale(y1_unit)

            ax1.plot(self.current_dist[0]*x_scale, self.current_dist[1]*y1_scale, ls='--',
                     lw=2, color='indigo')
            ax1.set_ylabel("$I$ " + y1_unit_label, fontsize=LABEL_FONT_SIZE, labelpad=LABEL_PAD)
            ax1.tick_params(labelsize=TICK_FONT_SIZE)

            plt.title(r"$\sigma_t$ = %s " % float("%.2g" % (self.St*x_scale))
                      + x_unit_label.replace('(', '').replace(')', '')
                      + r", $\sigma_\delta$ = %s " % float("%.2g" % self.Sdelta)
                      + r", $Q$ = %s pC" % float("%.2g" % (self.charge*1e12)),
                      fontsize=TITLE_FONT_SIZE, y=1.02)

        if x_lim is not None:
            ax.set_xlim(x_lim)
        if y_lim is not None:
            ax.set_ylim(y_lim)

        plt.tight_layout()
        if output:
            image_file = os.path.join(os.path.dirname(self.particle_file), output)
            plt.savefig(image_file, dpi=dpi)
            print((image_file + ' generated.'))
        else:
            plt.show()

    @staticmethod
    def get_column(name):
        """Return the column name in PhaseSpace object"""
        if name == 'xp':
            return 'px'
        elif name == 'yp':
            return 'py'
        else:
            return name

    @staticmethod
    def sample_data(x, y, bins=None, sigma=None, sample=20000):
        """Sample the data and calculate the density map

        Parameters
        ----------
        x: pandas.Series
            x data
        y: pandas.Series
            y data
        bins: int or [int, int]
            No. of bins used in numpy.histogram2d().
        sigma: numeric
            Standard deviation of Gaussian kernel of the Gaussian filter.
        sample: scalar >=0
            If sample < 1.0, sample by fraction;
            else, sample by count (round to integer).

        Returns
        -------
        x_sample: pandas.Series
            sampled x data.
        y_sample: pandas.Series
            sampled y data
        z: numpy.ndarray.
            Normalized density at each sample point.
        """
        if int(sample) > 0:
            n = int(sample)
        elif sample > 0:
            n = int(sample*len(x))
        else:
            raise ValueError("Negative sample value!")

        H, x_edges, y_edges = np.histogram2d(x, y, bins=bins)
        x_center = (x_edges[1:] + x_edges[0:-1]) / 2
        y_center = (y_edges[1:] + y_edges[0:-1]) / 2
        H_blurred = gaussian_filter(H, sigma=sigma)

        if len(x) > n:
            i_sample = random.sample(list(range(len(x))), n)
            x_sample = x.iloc[i_sample]
            y_sample = y.iloc[i_sample]
        else:
            i_sample = np.array(list(range(len(x))))
            x_sample = x
            y_sample = y

        posx = np.digitize(x_sample, x_center)
        posy = np.digitize(y_sample, y_center)
        z = H_blurred[posx - 1, posy - 1]
        z = z/z.max()

        return x_sample, y_sample, z, i_sample


class LinePlot(BeamEvolution):
    """Inherit from BeamEvolution object"""
    def plot(self, names, x_unit=None, y_unit=None, x_lim=None, y_lim=None,
             lattice_file=None, figsize=None, colors=None, styles=None,
             output='', dpi=300):
        """Plot the beam parameter evolution along the beamline

        Parameters
        ----------
        names: string/list/tuple
            A variable or a list of variables to plot
        x_unit: string
            Units of y1 axis.
        y_unit: string
            Units of y1 axis.
        x_lim -> tuple
            Limits for the x axis as (left, right)
        y_lim -> tuple
            Limits for the y axis as (left, right)
        lattice_file: string
            File contains the magnetic lattice layout.
        figsize: None/tuple
            Size of the figure (width, height).
        colors: None/string/list/tuple
            A color or a list of colors.
        styles: None/string/list/tuple
            A line style or a list of line styles
        output: string
            Name of the output file.
        dpi: int
            dpi of the png output.
        """
        names = _to_tuple(names)

        if len(names) > 2:
            raise ValueError("\nToo many (>2) variables to plot at one axis!")

        if not x_unit:
            x_unit = 'm'

        line_data = self.data
        if line_data is None:
            print("\nError: NoneType data!")
            raise SystemExit

        x_unit_label, x_scale = unit_scale(x_unit)

        if figsize is None:
            figsize = (8, 5)
        fig, ax = plt.subplots(figsize=figsize)
        ax.margins(AX_MARGIN)
        ax.xaxis.set_major_locator(
            ticker.MaxNLocator(MAX_LOCATOR, symmetric=False))
        ax.yaxis.set_major_locator(
            ticker.MaxNLocator(MAX_LOCATOR, symmetric=False))

        if colors is None:
            colors = ['dodgerblue', 'firebrick']
        else:
            colors = _to_tuple(colors)

        if styles is None:
            styles = ['-', '--']
        else:
            styles = _to_tuple(styles)

        # Set the default y unit
        if y_unit is None:
            y_unit = default_unit(names[0])

        y_unit_label, y_scale = unit_scale(y_unit)

        name_labels = []
        for i in range(len(names)):
            name_labels.append(name2label(names[i]))
            try:
                ax.plot(line_data['z']*x_scale, line_data[names[i]]*y_scale,
                        c=colors[i], ls=styles[i], lw=2, label=name_labels[i])
            except KeyError as inst:
                print(("\nUnknown key: {}".format(inst)))
                print(("\nValid keys ares: {}".format(line_data.columns.values)))
                plt.close()  # Prevent showing failed plots afterwards
                raise SystemExit

        ax.set_xlabel("$z$ " + x_unit_label,
                      fontsize=LABEL_FONT_SIZE, labelpad=LABEL_PAD)
        ax.set_ylabel("$,$ ".join(name_labels) + " " + y_unit_label,
                      fontsize=LABEL_FONT_SIZE, labelpad=LABEL_PAD)
        ax.tick_params(labelsize=TICK_FONT_SIZE, pad=TICK_PAD)

        if x_lim is not None:
            ax.set_xlim(x_lim)
        if y_lim is not None:
            ax.set_ylim(y_lim)

        if lattice_file is not None:
            pass

        if len(names) > 1:
            plt.legend(loc=0, fontsize=LEGEND_FONT_SIZE)
        plt.tight_layout()

        if output:
            image_file = os.path.join(os.path.dirname(self.root_name), output)
            plt.savefig(image_file, dpi=dpi)
            print((image_file + ' generated.'))
        else:
            plt.show()


def _to_tuple(input):
    """Convert input to a tuple"""
    if isinstance(input, tuple):
        output = input
    elif isinstance(input, list):
        output = tuple(input)
    elif isinstance(input, str):
        output = tuple([input])
    else:
        raise TypeError()

    return output


def name2label(name):
    """Return the label of a variable used in the plot

    Parameters
    ----------
    name: string
        Variable name

    return
    ------
        label of the variable
    """
    if name == 'gamma':
        return r"$\gamma$"
    elif name == 'SdE':
        return r"$\sigma_E$"
    elif name == 'Sx':
        return "$\sigma_x$"
    elif name == 'Sy':
        return "$\sigma_y$"
    elif name == 'Sz':
        return "$\sigma_z$"
    elif name == 'betax':
        return r"$\beta_x$"
    elif name == 'betay':
        return r"$\beta_y$"
    elif name == 'alphax':
        return r"$\alpha_x$"
    elif name == 'alphay':
        return r"$\alpha_y$"
    elif name == 'emitx':
        return r"$\varepsilon_x$"
    elif name == 'emity':
        return r"$\varepsilon_y$"
    elif name == 'emitx_tr':
        return r"$\varepsilon_{x, trace}$"
    elif name == 'emity_tr':
        return r"$\varepsilon_{y, trace}$"
    elif name == 'xp':
        return r"$x^\prime$"
    elif name == 'yp':
        return r"$y^\prime$"
    else:
        return r"${}$".format(name)


def default_unit(name):
    """Return the default unit of an variable"""
    if name == 'x' or name == 'y':
        return 'mm'
    elif name == 'z':
        return 'm'
    elif name == 'xp' or name == 'yp':
        return 'mrad'
    elif name == 't':
        return 'fs'
    elif re.match('beta', name):
        return 'm'
    elif name == 'SdE':
        return 'keV'
    elif name == 'Sx' or name == 'Sy':
        return 'mm'
    elif name == 'Sz':
        return 'um'
    elif re.match('emit', name):
        return 'um'
    elif name == 'St':
        return 'fs'
    elif name == 'p':
        return 'mc'
    else:
        return ''


def unit_scale(unit):
    """Obtain the label and scaling factor of the unit

    Parameter
    ---------
    unit: string
        Name of the unit

    Returns
    -------
    unit_label: string:
        label of the unit
    scale: int/float
        Scaling factor of the unit
    """
    unit_label = unit

    if unit == 'MeV':
        scale = 1.0e-6
    elif unit in ['kA', 'keV']:
        scale = 1.0e-3
    elif unit == '' or unit in ['m', 'rad', 's', 'A', 'mc']:
        scale = 1.0
    elif unit in ['mm', 'mrad', 'ms']:
        scale = 1.0e3
    elif unit in ['um', 'urad', 'us']:
        scale = 1.0e6
        if unit == 'um':
            unit_label = '$\mu$m'
        elif unit == 'urad':
            unit_label = '$\mu$rad'
        elif unit == 'us':
            unit_label = '$\mu$s'
    elif unit == 'ps':
        scale = 1.0e12
    elif unit == 'fs':
        scale = 1.0e15
    else:
        raise ValueError("\nUnknown unit!")

    if unit_label:
        unit_label = "(" + unit_label + ")"

    return unit_label, scale


if __name__ == "__main__":
    # Test
    # psplot = PhaseSpacePlot('test/injector.0600.001', 'astra')
    psplot = PhaseSpacePlot('examples/plots/fort.140', 'impact', 3.1e-12, cut_tail=0.1)
    psplot.plot('t', 'x', density_plot=False, alpha=0.5)
    psplot.plot('t', 'p')
    psplot.plot('x', 'xp', x_unit='um', density_plot=False)
    psplot.plot('x', 'xp', x_unit='um', output='x-xp.png')

    # lineplot = LinePlot('test/injector', 'astra')
    # lineplot = LinePlot('examples/plots/fort', 'impact')
    # lineplot.plot('Sz')
    # lineplot.plot(['betax', 'betay'], output='betaxy.png')
    # lineplot.plot(['Sx', 'Sy'], colors=('red', 'blue'), styles=(':', '-.'))
    # lineplot.plot(['Sx', 'Sy'], x_lim=[30, 31], y_lim=(0, 0.25))
