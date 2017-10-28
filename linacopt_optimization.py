#!/usr/bin/python

"""

Author: Jun Zhu

"""
from linacopt_variable import Variable
from linacopt_covariable import CoVariable
from linacopt_staticvariable import StaticVariable
from linacopt_objective import Objective
from linacopt_constraint import Constraint

import pyOpt


class Optimization(pyOpt.Optimization):
    """Inherited from pyOpt.Optimization class"""
    def __init__(self, name, obj_fun, *args, **kwargs):
        """Initialize Optimization object"""
        self._covariables = {}
        self._staticvariables = {}

        super(Optimization, self).__init__(name, obj_fun, *args, **kwargs)

    def add_obj(self, *args, **kwargs):
        print("Deprecated! Please use set_obj method")
        self.set_obj(*args, **kwargs)

    def set_obj(self, *args, **kwargs):
        """Add a objective into the objective set"""
        i = self.firstavailableindex(self._objectives)

        if (len(args) > 0) and isinstance(args[0], Objective):
            self._objectives[i] = args[0]
        else:
            try:
                self._objectives[i] = Objective(*args, **kwargs)
            except:
                raise ValueError("Input is not a valid for a Objective "
                                 "Object instance\n")

        self._remove_duplicity(self._objectives, i)

    def del_obj(self, ix):
        """Delete the objective *ix* from the objective set

        arguments:
        - ix -> int/string: index (if int) or name (if string)
        """
        index = -1
        if isinstance(ix, int):
            if ix < 0:
                raise ValueError("Index must be an integer >= 0.\n")
            elif ix not in self._objectives:
                raise ValueError("Index is not found.\n")
            else:
                index = ix
        elif isinstance(ix, str):
            index = self._find_index_by_name(ix, self._objectives)
            if index < 0:
                raise ValueError("Objective '{}' is not found.\n".format(ix))

        del self._objectives[index]

        self.shift_index(self._objectives, index)

    def get_objset(self):
        """Get the objective set"""
        return self._objectives

    def add_con(self, *args, **kwargs):
        print("Deprecated! Please use set_con method")
        self.set_con(*args, **kwargs)

    def set_con(self, *args, **kwargs):
        """Add a constraint into the constraint set"""
        i = self.firstavailableindex(self._constraints)

        if (len(args) > 0) and isinstance(args[0], Constraint):
            self._constraints[i] = args[0]
        else:
            try:
                self._constraints[i] = Constraint(*args, **kwargs)
            except:
                raise ValueError("Input is not a valid for a Constraint "
                                 "Object instance\n")

        self._remove_duplicity(self._constraints, i)

    def del_con(self, ix):
        """Delete the constraint *ix* from the constraint set

        arguments:
        - ix -> int/string: index (if int) or name (if string)
        """
        index = -1
        if isinstance(ix, int):
            if ix < 0:
                raise ValueError("Index must be an integer >= 0.\n")
            elif ix not in self._constraints:
                raise ValueError("Index is not found.\n")
            else:
                index = ix
        elif isinstance(ix, str):
            index = self._find_index_by_name(ix, self._constraints)
            if index < 0:
                raise ValueError("Constraint '{}' is not found.\n".format(ix))

        del self._constraints[index]

        self.shift_index(self._constraints, index)

    def get_conset(self):
        """Get the constraint set"""
        return self._constraints

    def add_var(self, *args, **kwargs):
        print("Deprecated! Please use set_var method")
        self.set_var(*args, **kwargs)

    def set_var(self, *args, **kwargs):
        """Add a variable into the variable set"""
        i = self.firstavailableindex(self._variables)

        if (len(args) > 0) and isinstance(args[0], Variable):
            self._variables[i] = args[0]
        else:
            try:
                self._variables[i] = Variable(*args, **kwargs)
            except:
                raise ValueError(
                    "Input is not valid for a Variable Object instance.\n")

        # Try to inherit the value from the co-variable set.
        ii = self._find_index_by_name(self._variables[i].name,
                                      self._covariables)
        if ii >= 0:
            self._variables[i].value = self._covariables[ii].value
            self.del_covar(ii)
            print "\n{} was changed from co-variable to variable!". \
                format(self._variables[i].name)

        # Try to inherit the value from the static variable set.
        ii = self._find_index_by_name(self._variables[i].name,
                                      self._staticvariables)
        if ii >= 0:
            self._variables[i].value = self._staticvariables[ii].value
            self.del_staticvar(ii)
            print "\n{} was changed from static variable to variable!".\
                format(self._variables[i].name)

        self._remove_duplicity(self._variables, i)

    def del_var(self, ix):
        """Delete the variable *ix* from the variable set

        arguments:
        - ix -> int/string: index (if int) or name (if string)
        """
        index = -1
        if isinstance(ix, int):
            if ix < 0:
                raise ValueError("Index must be an integer >= 0.\n")
            elif ix not in self._variables:
                raise ValueError("Index is not found.\n")
            else:
                index = ix
        elif isinstance(ix, str):
            index = self._find_index_by_name(ix, self._variables)
            if index < 0:
                raise ValueError("Variable '{}' is not found.\n".format(ix))

        del self._variables[index]

        self.shift_index(self._variables, index)

    def get_varset(self):
        """Get the variable set"""
        return self._variables

    def add_covar(self, *args, **kwargs):
        print("Deprecated! Please use set_covar method")
        self.set_covar(*args, **kwargs)

    def set_covar(self, *args, **kwargs):
        """Add a co-variable into the co-variable set"""
        i = self.firstavailableindex(self._covariables)

        if (len(args) > 0) and isinstance(args[0], CoVariable):
            self._covariables[i] = args[0]
        else:
            try:
                self._covariables[i] = CoVariable(*args, **kwargs)
            except:
                raise ValueError(
                    "Input is not valid for a CoVariable object instance.\n")

        # Variable and static-variable cannot be changed to co-variable since
        # they will lose their values in the last optimization.
        try:
            self.del_var(self._covariables[i].name)
            print "\n{} was changed from variable to co-variable!". \
                format(self._covariables[i].name)
            print "Warning: changing variable to co-variable may lose "\
                  "the optimized result!\n"
        except ValueError:
            pass

        try:
            self.del_staticvar(self._covariables[i].name)
            print "\n{} was changed from static variable to co-variable!". \
                format(self._covariables[i].name)
            print "Warning: changing static variable to co-variable may lose "\
                  "the optimized result!\n"
        except ValueError:
            pass

        self._remove_duplicity(self._covariables, i)

    def del_covar(self, ix):
        """Delete the co-variable *ix* from the co-variable set

        arguments:
        - ix -> int/string: index (if int) or name (if string)
        """
        index = -1
        if isinstance(ix, int):
            if ix < 0:
                raise ValueError("Index must be an integer >= 0.\n")
            elif ix not in self._covariables:
                raise ValueError("Index is not found.\n")
            else:
                index = ix
        elif isinstance(ix, str):
            index = self._find_index_by_name(ix, self._covariables)
            if index < 0:
                raise ValueError("Co-variable '{}' is not found.\n".format(ix))

        del self._covariables[index]

        self.shift_index(self._covariables, index)

    def get_covarset(self):
        """Get co-variable set"""
        return self._covariables

    def add_staticvar(self, *args, **kwargs):
        print("Deprecated! Please use set_staticvar method")
        self.set_staticvar(*args, **kwargs)

    def set_staticvar(self, *args, **kwargs):
        """Add a static variable into the static variable set"""
        i = self.firstavailableindex(self._staticvariables)

        if (len(args) > 0) and isinstance(args[0], StaticVariable):
            self._staticvariables[i] = args[0]
        else:
            try:
                self._staticvariables[i] = StaticVariable(*args, **kwargs)
            except:
                raise ValueError(
                    "Input is not valid for a StaticVariable object instance\n")

        # Try to inherit the value from the variable set.
        ii = self._find_index_by_name(self._staticvariables[i].name,
                                      self._variables)
        if ii >= 0:
            self._staticvariables[i].value = self._variables[ii].value
            self.del_var(ii)
            print "\n{} was changed from variable to static variable!\n". \
                format(self._staticvariables[i].name)

        # Try to inherit the value from the co-variable set.
        ii = self._find_index_by_name(self._staticvariables[i].name,
                                      self._covariables)
        if ii >= 0:
            self._staticvariables[i].value = self._covariables[ii].value
            self.del_covar(ii)
            print "\n{} was changed from co-variable to static variable!\n". \
                format(self._staticvariables[i].name)

        if self._staticvariables[i].value is None:
            raise ValueError("Unknown value of static variable: {}".
                             format(self._staticvariables[i].name))

        self._remove_duplicity(self._staticvariables, i)

    def del_staticvar(self, ix):
        """Delete static variable *ix* from the static variable set

        arguments:
        - ix -> int/string: index (if int) or name (if string)
        """
        index = -1
        if isinstance(ix, int):
            if ix < 0:
                raise ValueError("Index must be an integer >= 0.\n")
            elif ix not in self._staticvariables:
                raise ValueError("Index is not found.\n")
            else:
                index = ix
        elif isinstance(ix, str):
            index = self._find_index_by_name(ix, self._staticvariables)
            if index < 0:
                raise ValueError("Static-variable '{}' is not found.\n".
                                 format(ix))

        del self._staticvariables[index]

        self.shift_index(self._staticvariables, index)

    def get_staticvarset(self):
        """Get the static variable set"""
        return self._staticvariables

    @staticmethod
    def shift_index(dict_, index):
        """"""
        while len(dict_) > index:
            dict_[index] = dict_[index+1]
            del dict_[index+1]
            index += 1

    @staticmethod
    def _remove_duplicity(dict_, ii):
        """Remove duplicity

        If the name of the ii-th item is equal to the name of another
        item with index i, the i-th item is set equal to the ii-th item
        and the latter is then removed.
        """
        for i, item in dict_.iteritems():
            if i != ii and item.name == dict_[ii].name:
                dict_[i] = dict_[ii]
                del dict_[ii]
                return

    @staticmethod
    def _find_index_by_name(name, dict_):
        """Return the index of object with object.name == name"""
        index = -1
        for key, ele in dict_.iteritems():
            if name == ele.name:
                index = key

        return index

    def __str__(self):
        """Print Structured Optimization Problem

        Overwrite the original method to include the summary of
        co-variables, static variables as well as functions in
        constraints and objectives.
        """
        text = "\nOptimization Problem -- %s\n%s\n" % (self.name, '='*80)

        text += "\nObjectives:\n"\
                "  {:18}  {:11}  {:11}  {:16}\n"\
                .format('Name', 'Value', 'Optimum', 'Function')

        for index in self._objectives.keys():
            lines = str(self._objectives[index]).split('\n')
            text += lines[1] + '\n'

        if len(self._constraints.keys()) > 0:
            text += "\nConstraints:\n"\
                    "  {:18}  {:11}  {:48}\n"\
                    .format('Name', 'Value', 'Bound')
            for index in self._constraints.keys():
                lines = str(self._constraints[index]).split('\n')
                text += lines[1] + '\n'

        text += "\nVariables (c - continuous, i - integer, d - discrete):\n"\
                "  {:18}  {:6}  {:11}  {:11}  {:11}\n"\
                .format('Name', 'Type', 'Value', 'Lower Bound', 'Upper Bound')
        for index in self._variables.keys():
            lines = str(self._variables[index]).split('\n')
            text += lines[1] + '\n'

        if len(self._covariables.keys()) > 0:
            text += "\nCo-variables:\n"\
                    "  {:18}  {:16}  {:11}  {:11}\n"\
                    .format('Name', 'Dependent(s)', 'Slope(s)', 'Intercept')
            for index in self._covariables.keys():
                lines = str(self._covariables[index]).split('\n')
                text += lines[1] + '\n'

        if len(self._staticvariables.keys()) > 0:
            text += "\nStatic variables:\n"\
                    "  {:18}  {:11}\n".format('Name', 'Value')
            for index in self._staticvariables.keys():
                lines = str(self._staticvariables[index]).split('\n')
                text += lines[1] + '\n'

        return text
