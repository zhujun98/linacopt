#!/usr/bin/python

"""
linacopt_staticvariable

Hold the StaticVariable class

Tested on:
__________
Ubuntu 14.04

Developer:
__________
- Jun Zhu

History:
________
Last modified on 30/12/2016

"""


class StaticVariable(object):
    """Optimization StaticVariable class

    Static variable is useful when running concatenate simulations.
    """

    def __init__(self, name, value=None):
        """Initialize StaticVariable object

        Parameters
        ----------
        name: string
            Name of the static variable.
        value: int/float
            Value of the static variable. The value will be ignored
            if a variable with the same name is found in the variable
            set or co-variable set of the solution.
        """
        self.name = name
        self.value = value

    def __str__(self):
        """Print structured list of static variables.

        Overwrite the original method.
        """
        return ('  {:18}  {:11}\n'.format('Name', 'Value')
                + '  {:18}  {:11.4e}\n'.format(self.name, self.value))

