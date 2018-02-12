#!/usr/bin/python

"""

Author: Jun Zhu

"""
import os

from .data_processing import PhaseSpace


class FitPoints(object):
    """Object with attributes being PhaseSpace objects"""
    def __init__(self, particle_type):
        self.particle_type = particle_type

    def set_point(self, name, particle_file, particle_type=None, **kwargs):
        """Add a PhaseSpace object as an attribute

        :param name: string
            Name of the new attribute.
        :param particle_file: string
            Name of the particle file.
        :param particle_type: string
            Type of the particle file.

        Additional keyword arguments are passed to PhaseSpace.__init__()
        """
        try:
            super().__delattr__(name)
            print("\nWarning: point {} will be replaced!".format(name))
        except AttributeError:
            pass

        if particle_type is None:
            particle_type = self.particle_type

        this_point = PhaseSpace(particle_file, particle_type, opt=True, **kwargs)
        super().__setattr__(name, this_point)

    def del_point(self, name):
        """Delete a point by name.

        :param name: string
            Name of the new attribute.
        """
        super().__delattr__(name)

    def update(self, **kwargs):
        """Update all the attributes which are PhaseSpace objects"""
        for value in list(self.__dict__.values()):
            if isinstance(value, PhaseSpace):
                value.update(**kwargs)

    def __str__(self):
        """Print structured list of attributes"""
        text = '\nFit points:\n'

        text += '  {:18}  {:16}  {:12}  {:12}  {:12}\n'.format(
            'Name', 'File', 'Type', 'Cut halo', 'Cut tail')
        for key, value in self.__dict__.items():
            if isinstance(value, PhaseSpace):
                text += '  {:18}  {:16}  {:12}  {:12}  {:12}\n'.format(
                    key, os.path.basename(value.particle_file), value.particle_type,
                    str(value.cut_halo), str(value.cut_tail))

        return text
