#!/usr/bin/python

"""

Author: Jun Zhu

"""
from abc import ABC
import re


class LinacOptData(ABC):
    """Data object used in linacopt optimizer"""
    def __init__(self, particle_type=None):
        """"""
        self._particle_type = None
        if particle_type is not None:
            self.particle_type = particle_type
        else:
            self.particle_type = 'astra'

    @property
    def particle_type(self):
        """Indicates the code that generated the data"""
        return self._particle_type

    @particle_type.setter
    def particle_type(self, value):
        """Setter"""
        if re.search(r'^astra', value, re.IGNORECASE):
            self._particle_type = 'astra'
        elif re.search(r'^impact', value, re.IGNORECASE):
            self._particle_type = 'impact'
        else:
            raise ValueError(
                "The particle_type must be 'astra' or 'impact'!\n")
