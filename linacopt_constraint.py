#!/usr/bin/python

"""
Note: the lower and upper keyword arguments in the original Constraint
class do not work!


Author: Jun Zhu

"""
from pyOpt import Constraint as pyOptConstraint


INF = 10.E+20


class Constraint(pyOptConstraint):
    """Inherited from pyOpt.Constraint class.

    One can specify one of the three constraint formats:
    [-INF, upper], [lower, INF] and [equal - tol, equal + tol]
    """
    def __init__(self, name, func=None, **kwargs):
        """Initialize Constraint object.

        Parameters
        ----------
        func: function object
            Name of the constraint function.
        upper: float
            Upper boundary, ignored if lower if specified.
        lower: float
            Lower boundary.
        equal: float
            Equal boundary.
        tol: float
            Tolerance for equal boundary; negative number for percent,
            positive number for absolute value.

        Note: the attributes 'lower' and 'upper' in the pyOpt.constraint
        object does not take effect. In the new Constraint object,
        the normalize() method will convert the value corresponding to
        a constraint type of (-INF, 0] which can then pass to pyOpt.
        """
        self.name = name
        self.value = 0.0
        # When pyOpt see the object, it will always consider this is an
        # inequality constraint between -inf and 0.
        self.type = 'i'
        self.upper = 0.0
        self.lower = -1.0*float(INF)

        if not callable(func):
            raise TypeError('{} is not a function!\n'.format(func))
        else:
            self.func = func

        self._upper_ = None
        self._lower_ = None
        self._equal_ = None
        self._tol = 0.0
        self._setup_boundary(**kwargs)

    @property
    def upper_(self):
        return self._upper_

    @upper_.setter
    def upper_(self, value):
        self._upper_ = value
        self._lower_ = -float(INF)

    @property
    def lower_(self):
        return self._lower_

    @lower_.setter
    def lower_(self, value):
        self._lower_ = value
        self._upper_ = float(INF)

    @property
    def equal_(self):
        return self._equal_

    @equal_.setter
    def equal_(self, value):
        self._equal_ = value

    def _setup_boundary(self, **kwargs):
        """Setup the boundaries of the constraint."""
        for key in kwargs:
            if key == 'upper':
                self.upper_ = kwargs[key]
                break

            if key == 'lower':
                self.lower_ = kwargs[key]
                break

            if key == 'equal':
                self.equal_ = kwargs[key]
                break

        # Default value
        if [self.upper_, self.lower_, self.equal_] == [None]*3:
            self.upper_ = 0.0

        if self.equal_ is not None:
            for key in kwargs:
                if key == 'tol':
                    self._tol = kwargs[key]

            if self._tol < 0:
                self._tol = -1.0*self._tol*self.equal_
            else:
                pass

    def update(self, *args, **kwargs):
        """Update the value attribute"""
        self.value = self.func(*args, **kwargs)

    def normalize(self):
        """Calculate and return the value used by pyOpt"""
        if self.equal_ is not None:
            g_normalized = abs(self.value - self.equal_) - self._tol
        else:
            # self.lower_ == -INF
            if abs(self.lower_ + INF) < 1e-6:
                g_normalized = self.value - self.upper_
            else:
                g_normalized = self.lower_ - self.value

        return g_normalized

    def __str__(self):
        """Print structured list of constraints.

        Overwrite the original method.
        """
        if self.equal_ is None:
            bound = "{:.4e} <= {}() <= {:.4e}".format(
                self.lower_, self.func.__name__, self.upper_)
        else:
            bound = "{1:.4e} - {2:.4e} <= {0}() <= {1:.4e} + {2:.4e}".\
                format(self.func.__name__, self.equal_, self._tol)

        return ('  {:18}  {:11}  {:48}\n'.
                format('Name', 'Value', 'Bound')
                + '  {:18}  {:11.4e}  {:48}\n'.
                format(self.name, self.value, bound))
