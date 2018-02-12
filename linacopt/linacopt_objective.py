#!/usr/bin/python

"""

Author: Jun Zhu

"""
from pyOpt import Objective as pyOptObjective


class Objective(pyOptObjective):
    """Inherited from pyOpt.Objective class.

    Each Objective object has a function attribute which calculates
    its value.
    """
    def __init__(self, name, func=None, optimum=-1.0e21):
        """Initialize Objective object

        :param func: function object
            Objective function.
        :param optimum: None/float
            The object of the objective.
        """
        if not callable(func):
            raise TypeError('{} is not a function!\n'.format(func))
        else:
            self.func = func

        self.name = name
        self.value = 0.0

        self.optimum = optimum

    def update(self, *args, **kwargs):
        """Update the value attribute."""
        self.value = self.func(*args, **kwargs)

    def normalize(self):
        """Return the objective seen by pyOpt"""
        if self.value > self.optimum:
            return self.value
        else:
            return self.optimum

    def __str__(self):
        """Print structured list of objectives.

        Overwrite the original method.
        """
        return ('  {:18}  {:12}  {:12}  {:16}\n'.
                format('Name', 'Value', 'optimum', 'Function')
                + '  {:18}  {:11.4e}  {:11.4e}  {:16}\n'.
                format(self.name, self.value, self.optimum, self.func.__name__))
