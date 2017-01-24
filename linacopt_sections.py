#!/usr/bin/python

"""
linacopt_sections

Hold one class.

    - Sections: hold attributes being BeamEvolution objects

Tested on:
__________
Ubuntu 14.04
Ubuntu 16.04

Developer:
__________
- Jun Zhu

History:
________
Last modified on 22/01/2017

"""
import os

from beam_evolutions import BeamEvolution


class Sections(object):
    """Object with attributes being BeamEvolution objects"""
    def __init__(self, particle_type, root_name):
        self.particle_type = particle_type
        self.root_name = root_name

    def add_section(self, *args, **kwargs):
        print("Deprecated! Please use set_section method")
        self.set_section(*args, **kwargs)

    def set_section(self, name, root_name=None, particle_type=None, **kwargs):
        """Add a BeamEvolution object as an attribute

        Parameters
        ----------
        name: string
            Name of the new attribute.
        root_name: string
            The root name of the output files.
        particle_type: string
            Type of the particle file.

        Additional keyword arguments are passed to BeamEvolution.__init__().
        """
        try:
            self.__delattr__(name)
            print "\nWarning: section {} will be replaced!".format(name)
        except AttributeError:
            pass

        if root_name is None:
            root_name = self.root_name

        if particle_type is None:
            particle_type = self.particle_type

        this_section = BeamEvolution(root_name, particle_type, opt=True, **kwargs)
        super(Sections, self).__setattr__(name, this_section)

    def del_section(self, name):
        """Delete a section by name.

        Parameters
        ----------
        name: string
            Name of the new attribute.
        """
        super(Sections, self).__delattr__(name)

    def update(self):
        """Update all the attributes which are BeamEvolution objects"""
        for value in self.__dict__.values():
            if isinstance(value, BeamEvolution):
                value.update()

    def __str__(self):
        """Print structured list of attributes"""
        text = '\nSections:\n'

        text += '  {:18}  {:16}  {:12}  {:16}\n'.format(
            'Name', 'Root name', 'Type', 'z_lim (m)')

        for key, value in self.__dict__.iteritems():
            if isinstance(value, BeamEvolution):
                text += '  {:18}  {:16}  {:12}  {:16}\n'.format(
                    key, os.path.basename(value.root_name), value.particle_type, value.z_lim)

        return text