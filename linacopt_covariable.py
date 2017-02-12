#!/usr/bin/python

"""
linacopt_covariable

Hold the CoVariable class.

Tested on:
__________
Ubuntu 14.04
Ubuntu 16.04

Developer:
__________
- Jun Zhu

History:
________
Last modified on 12/01/2017

"""


class CoVariable(object):
    """Optimization CoVariable class."""
    def __init__(self, name, var_dp, slope=1.0, intercept=0.0):
        """Initialize CoVariable object

        The value of the variable is calculated by:
        covar = [slope].*[value of var_dp] + intercept

        Parameters
        ----------
        name: string
            Name of the co-variable.
        var_dp: string/list
            Name(s) of the dependent variable(s).
        slope: int/float/list
            Coefficient(s) in calculation.
        intercept: int/float
            Coefficient in calculation.
        """
        self.name = name

        if isinstance(var_dp, list) or isinstance(var_dp, tuple):
            self.var_dp = var_dp
        else:
            self.var_dp = [var_dp]

        self.intercept = intercept

        if isinstance(slope, list) or isinstance(slope, tuple):
            self.slope = slope
        else:
            self.slope = [slope]

        self.value = None

        if name in var_dp:
            raise ValueError(
                "The co-variable and dependent variable are the same.\n")

        if len(self.slope) != len(self.var_dp):
            raise ValueError(
                "slope does not have the same length as var_dp!\n")

    def __str__(self):
        """Print structured list of co-variables.

        Overwrite the original method.
        """
        if len(self.var_dp) == 1:
            return ('  {:18}  {:16} {:11}  {:11}\n'.
                    format('Name', 'Dependent', 'Slope', 'Intercept')
                    + '  {:18}  {:16}  {:11.4e}  {:11.4e}\n'.
                    format(self.name, self.var_dp[0], self.slope[0], self.intercept))
        else:
            return ('  {:18}  {:16}  {:11}  {:11}\n'.
                    format('Name', 'Dependents', 'Slopes', 'Intercept')
                    + '  {:18}  {:16}... {:11.4e}... {:11.4e}\n'.
                    format(self.name, self.var_dp[0], self.slope[0], self.intercept))
