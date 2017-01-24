#!/usr/bin/python

"""
linacopt_covariable

Hold the CoVariable class.

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


class CoVariable(object):
    """Optimization CoVariable class."""
    def __init__(self, name, var_dp, **kwargs):
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
        self.value = None

        if isinstance(var_dp, list) or isinstance(var_dp, tuple):
            self.var_dp = var_dp
        else:
            self.var_dp = [var_dp]

        if name in var_dp:
            raise ValueError(
                "The co-variable and dependent variable are the same.\n")

        self.slope_ = [1.0]*len(var_dp)
        self.intercept_ = 0.0
        for key in kwargs:
            if key == 'slope':
                if isinstance(kwargs[key], list) or \
                        isinstance(kwargs[key], tuple):
                    self.slope_ = kwargs[key]
                else:
                    self.slope_ = [kwargs[key]]

                if len(self.slope_) != len(self.var_dp):
                    raise ValueError(
                        "slope does not have the same length as var_dp!\n")

            elif key == 'intercept':
                self.intercept_ = kwargs[key]

    def __str__(self):
        """Print structured list of co-variables.

        Overwrite the original method.
        """
        if len(self.var_dp) == 1:
            return ('  {:18}  {:16} {:11}  {:11}\n'.
                    format('Name', 'Dependent', 'Slope', 'Intercept')
                    + '  {:18}  {:16}  {:11.4e}  {:11.4e}\n'.
                    format(self.name, self.var_dp[0], self.slope_[0], self.intercept_))
        else:
            return ('  {:18}  {:16}  {:11}  {:11}\n'.
                    format('Name', 'Dependents', 'Slopes', 'Intercept')
                    + '  {:18}  {:16}... {:11.4e}... {:11.4e}\n'.
                    format(self.name, self.var_dp[0], self.slope_[0], self.intercept_))
