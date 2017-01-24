#!/usr/bin/python

"""
linacopt_data.py
-------------

Hold one class
- LinacOptData
    Base class for objects storing data.

Developer:
__________
Jun Zhu

Tested on:
__________
Ubuntu 14.04

History:
________

Last modified on 30/12/2016

"""

import re


class LinacOptData(object):
    """Data object used in linacopt optimizer"""
    def __init__(self, particle_type):
        """"""
        self._particle_type = None
        if particle_type is not None:
            self.particle_type = particle_type
        else:
            raise ValueError("particle_type is NoneType!")

    @property
    def particle_type(self):
        """Indicates the code that generated the data"""
        return self._particle_type

    @particle_type.setter
    def particle_type(self, value):
        if re.search(r'^astra', value, re.IGNORECASE):
            self._particle_type = 'astra'
        elif re.search(r'^impact', value, re.IGNORECASE):
            self._particle_type = 'impact'
        else:
            raise ValueError(
                "The particle_type must be 'astra' or 'impact'!\n")