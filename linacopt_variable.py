#!/usr/bin/python

"""

Author: Jun Zhu

"""
from pyOpt import Variable as pyOptVariable


class Variable(pyOptVariable):
    """Inherited from pyOpt.Variable class"""
    def __init__(self, name, **kwargs):
        """Initialize Constraint object.

        Parameters
        ----------
        type: string
            Variable Type ('c'-continuous, 'i'-integer, 'd'-discrete),
            *Default* = 'c'
        value: int/float
            Variable Value, *Default* = 0.0
        upper: int/float
            Upper boundary.
        lower: int/float
            Lower boundary.
        """
        super(Variable, self).__init__(name, **kwargs)

    def __str__(self):
        """Print structured list of variables

        Overwrite the original method
        """
        return ('  {:18}  {:6}  {:11}  {:11}  {:11}\n'.
                format('Name', 'Type', 'Value', 'Lower Bound', 'Upper Bound')
                + '  {:18}  {:6}  {:11.4e}  {:11.4e}  {:11.4e}\n'.
                format(self.name, self.type, self.value, self.lower, self.upper))
