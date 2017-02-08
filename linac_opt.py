#!/usr/bin/python

"""
linac_Opt.py - A PYTHON script for optimizing linac.

Optimizers (SDPEN, ALPSO, NSGA2) in pyOpt are used in this
script to solve general constrained nonlinear optimization problems:

    min f(x) w.r.t. x

    s.t. g_j(x) = 0, j = 1, ..., m_e

        g_j(x) <= 0, j = m_e + 1, ..., m

        x_i_L <= x_i <= x_i_U, i = 1, ..., n

    where:

        x is the vector of design variables;
        f(x) is a nonlinear function;
        g(x) is a linear or nonlinear function;
        n is the number of design variables;
        m_e is the number of equality constraints;
        m is the total number of constraints (number of equality
        constraints: m_i = m - m_e).

Tested on:
__________
Ubuntu 14.04
Ubuntu 16.04

pyOpt version 1.2.0

Developer:
__________
- Jun Zhu

History:
________

version 1.00
Last modified on 24/01/2017
"""
import os
import subprocess
import re
import pickle
from glob import glob
from datetime import datetime

import pyOpt

from linacopt_optimization import Optimization
from linacopt_fitpoints import FitPoints
from linacopt_sections import Sections
from linacopt_data import LinacOptData
from beam_parameters import PhaseSpace
from beam_evolutions import BeamEvolution


INF = 1.0e21


class LinacOpt(LinacOptData):
    """LinacOpt class.

    Attributes
    ----------
    input_file: string
        Name of the input file.
    _input_template_lines: string
        Strings read from the input file template.
    particle_type: string
        Type of the particle file (corresponding code name).
    _prob_name: string
        Root name of log_file and solution_file.
    name: string
        Name of the optimization problem.
    log_file: string
        File name recording the optimization history of the current run.
    solution_file: string
        File name recording the solution of the current run.
    run_code: string
        String that can run the code in the shell.
    _time_consumption: float
        Time consumption (in second) of the last iteration.
    _n_iter: int
        Number of iterations in the current run.
    _n_restart: int
        Number of restart runs.
    _n_fail: int
        Number of consecutive fails.
    opt_prob: object
        Optimization object.
    optimizer: object
        Optimizer object in pyOpt.
    fit_points: FitPoints object
        Object with attributes being PhaseSpace objects.
    sections: Sections object
        Object with attributes being BeamEvolution objects.
    _f: array-like
        Objects.
    _f2pyopt: array-like
        Objective values seen by pyOpt.
    _x: array-like
        Variables.
    _xc: array-like
        Co-variables.
    _g2pyopt: array-like
        Constraint values seen by pyOpt.
    _g: array-like
        Constraint values calculated by func.
    _restart: None or int
        None for new run and int for restarting from the restart-th run.
    _run_once: Boolean
        True for stopping after one simulation. For debug.
    run_code: string
        String that can run the code in the shell.
    time_out: float/int
        Maximum allowed run time (in s) of one simulation.
    complete_shell: Boolean
        If True, the shell will only execute the string 'run_code',
        which means the input file and the time_out will be ignored.
    """
    def __init__(self, path_name=None, input_file=None, input_template=None,
                 particle_type=None, prob_name='opt_prob',
                 restart=None, run_once=False, max_fail=20):
        """Initialization.

        Parameters
        ----------
        path_name: string
            Path of the simulation directory.
        input_file: string
            Name of the input file.
        input_template: string
            Name of the input template file. The patterns (<string>) in
            the template file will be replaced at the beginning of
            each simulation to generate a new input file.
        particle_type: string
            Type of the particle file (corresponding code name).
        prob_name: string
            Root name of log_file and solution_file, default = 'opt_prob'.
        restart: None/int
            None for new run and int for restarting from the
            restart-th run.
        run_once: Boolean
            True for stopping after one simulation. For debug.
        max_fail: int
            Max number of allowed consecutive fails.
        """
        os.chdir(path_name)

        self.input_file = os.path.join(os.getcwd(), input_file)
        # input_template file will only be read once.
        self._input_template_lines = open(input_template).readlines()

        super(LinacOpt, self).__init__(particle_type)

        self._prob_name = prob_name
        self.name = prob_name + '-' + datetime.now().strftime('%Y-%m-%d-%H-%M-%S')
        self.log_file = None
        self.solution_file = None

        self.run_code = None
        self.time_out = 1200
        self.complete_shell = False

        self._time_consumption = 0.0
        self._n_iter = 0
        self._n_restart = 0
        self._n_fail = 0
        self._max_fail = max_fail

        # Set up the problem
        self.opt_prob = Optimization(self.name, self.obj_func)
        # Set up the optimizer
        self._optimizer_name = None
        self._optimizer = None

        if self.particle_type == 'impact':
            root_name = 'fort'
        else:
            root_name = os.path.basename(self.input_file).split('.')[0]

        self.fit_points = FitPoints(self.particle_type)
        self.sections = Sections(self.particle_type, root_name)

        # The following attributes(_f, _x, _g, _g2pyopt, _xc, _p) will
        # be updated after each iteration by the self.obj_func()
        self._f = None  # Object
        self._f2pyopt = None  # Object values seen by pyOpt
        self._x = None  # Variables
        self._g2pyopt = None  # Constraint values seen by pyOpt
        self._g = None  # Constraint values calculated by func
        self._xc = None  # Co-variables

        self._restart = restart
        if self._restart is not None:
            assert type(restart) == int and restart > 0

        self._run_once = run_once  # For debug

    @property
    def optimizer(self):
        """Property optimizer getter."""
        return self._optimizer

    @optimizer.setter
    def optimizer(self, name):
        """Property optimizer setter.

        Parameters
        ----------
        name: string
            Name of the optimizer.
        """
        if re.search(r'^sdpen', name, re.IGNORECASE):
            self._optimizer = pyOpt.SDPEN()
            self._optimizer_name = 'sdpen'
            self._optimizer.setOption('ifile', 'SDPEN.out.' +
                                      "{:03d}".format(self._n_restart))

        elif re.search(r'^alpso', name, re.IGNORECASE):
            self._optimizer = pyOpt.ALPSO()
            self._optimizer_name = 'alpso'

        elif re.search(r'^nsga', name, re.IGNORECASE):
            self._optimizer = pyOpt.NSGA2()
            self._optimizer_name = 'nsga2'

        else:
            raise ValueError("Unknown optimizer!\n")

    def set_optimizer(self, name):
        """Set the optimizer.

        Parameters
        ----------
        name: string
            Name of the optimizer.
        """
        self.optimizer = name

    def solve(self, run_code=None, time_out=None, complete_shell=None):
        """Run the optimization and print the result.

        Parameters
        ----------
        run_code: None/string
            String that can run the code in the shell.
        time_out: None/float/int
            Maximum allowed run time (in s) of one simulation.
        complete_shell: None/Boolean
            If True, the shell will only execute the string 'run_code',
            which means the input file and the time_out will be ignored.
        """
        # Continuous run
        # Read the optimized variable set from the .pkl file
        if self._restart is not None and self._restart > self._n_restart:
            pickle_file = self._prob_name + '.sol.{:03d}.pkl'.format(self._n_restart)
            with open(pickle_file, 'rb') as fp:
                last_sol = pickle.load(fp)
            self.name = last_sol[0]
            self.opt_prob.name = self.name
            self.opt_prob._variables = last_sol[1]
            self.opt_prob._covariables = last_sol[2]
            self.opt_prob._staticvariables = last_sol[3]
            self.opt_prob._objectives = last_sol[4]
            self.opt_prob._constraints = last_sol[5]
            self.run_code = last_sol[6]
            self.time_out = last_sol[7]
            self.complete_shell = last_sol[8]

            print("\n" + "*"*80 + "\n" + "Read the solution set from {}\n".
                  format(pickle_file) + "*"*80 + "\n")

            self._n_restart += 1
            return

        elif self._n_restart == 0:
            existing_log_files = glob(self._prob_name + '.log.*')
            existing_sol_files = glob(self._prob_name + '.sol.*')
            for this_file in existing_log_files:
                os.remove(this_file)
                print("File removed: {}".format(this_file))
            for this_file in existing_sol_files:
                os.remove(this_file)
                print("File removed: {}".format(this_file))

        if run_code is not None:
            self.run_code = run_code
        if time_out is not None:
            self.time_out = time_out
        if complete_shell is not None:
            self.complete_shell = complete_shell

        self.log_file = self._prob_name + '.log.{:03d}'.format(self._n_restart)
        self.solution_file = self._prob_name + '.sol.{:03d}'.format(self._n_restart)

        print "\n{}".format('*'*80)
        print "Start solving the following problem with " + \
              "pyOpt.{} on \n{} {}".format(self._optimizer_name.upper(),
                                           self.run_code, self.input_file)
        print "{}".format('*'*80)

        print self.opt_prob, self.fit_points, self.sections

        t0 = datetime.now()

        self.optimizer(self.opt_prob)

        self._run_optimized()

        self._time_consumption = datetime.now() - t0

        print self._solution_text()

        with open(self.solution_file, 'wb') as fp:
            fp.write(self._solution_text())

        # Store objects for continuous runs
        # ---------------------------------
        last_sol = [self.name,
                    self.opt_prob.get_varset(),
                    self.opt_prob.get_covarset(),
                    self.opt_prob.get_staticvarset(),
                    self.opt_prob.get_objset(),
                    self.opt_prob.get_conset(),
                    self.run_code,
                    self.time_out,
                    self.complete_shell]
        with open(self._prob_name + '.sol.{:03d}.pkl'.format(self._n_restart), 'wb') as fp:
            pickle.dump(last_sol, fp)

        self._n_iter = 0
        self._n_restart += 1

    def _run_optimized(self):
        """Run the simulation with the optimized variables.

        Update the values of the attributes in the Optimization object.
        """
        last_sol_key = max(self.opt_prob.getSolSet().keys())
        for key, item in self.opt_prob.getSolSet()[last_sol_key].getVarSet().iteritems():
            self.opt_prob.get_varset()[key].value = item.value
            self._x[key] = item.value

        self._update_covar()
        for key in self.opt_prob.get_covarset().keys():
            self.opt_prob.get_covarset()[key].value = self._xc[key]

        self._update_input()
        self._run_simulation()

        self._update_output()
        for key in self.opt_prob.get_objset().keys():
            self.opt_prob.get_objset()[key].value = self._f[key]
        for key in self.opt_prob.get_conset().keys():
            self.opt_prob.get_conset()[key].value = self._g[key]

        # Compare the result of pyOpt and the last simulation.
        tol = 1.0e-6
        for key in self.opt_prob.getSolSet()[last_sol_key].getObjSet().keys():
            old = self.opt_prob.getSolSet()[last_sol_key].getObjSet()[key].value
            new = self._f2pyopt[key]
            assert abs(new - old) < tol*max(abs(new), abs(old)), \
                "Different old and new results of objective functions!"

    def obj_func(self, x):
        """Objective function.

        Parameters
        ----------
        x: array-like
            Variables updated after each iteration.
        """
        self._n_iter += 1

        t0 = datetime.now()

        self._x = list(x)  # Use different address!!!

        self._update_covar()
        self._update_input()

        # t0_code = datetime.now()
        self._run_simulation()
        # dt_code = (datetime.now() - t0_code).total_seconds()
        # print dt_code

        nobj = len(self.opt_prob.get_objset().values())
        ncon = len(self.opt_prob.get_conset().values())
        try:
            self._update_output()
            self._n_fail = 0
            fail = 0

        except (IOError, ValueError) as e:
            print str(e)
            self._f = [INF]*nobj
            self._f2pyopt = [INF]*nobj
            self._g = [INF]*ncon
            self._g2pyopt = [INF]*ncon
            self._n_fail += 1
            fail = 1

        if self._n_fail > self._max_fail:
            raise IOError("Number of consecutive fails is larger than {}! "
                          "Please check the input file!\n".format(self._max_fail))

        dt = (datetime.now() - t0).total_seconds()

        self._write_history(dt)

        if self._run_once is True:
            raise SystemExit("Debugging: stopped after one pass!")

        return self._f2pyopt[0], self._g2pyopt, fail

    def _update_output(self):
        """Update the attributes used for optimization.

        update _f, _g, _g2pyopt, _p
        """
        self.fit_points.update()
        self.sections.update()

        self._f = []
        self._f2pyopt = []
        for obj in self.opt_prob.get_objset().values():
            try:
                obj.update(self.fit_points, self.sections)
            except TypeError:
                obj.update(self.fit_points)

            self._f.append(obj.value)
            self._f2pyopt.append(obj.normalize())

        self._g = []
        self._g2pyopt = []
        for con in self.opt_prob.get_conset().values():
            try:
                con.update(self.fit_points, self.sections)
            except TypeError:
                con.update(self.fit_points)

            self._g.append(con.value)
            self._g2pyopt.append(con.normalize())

    def _update_covar(self):
        """Update the co-variables."""
        xc = []
        for covar in self.opt_prob.get_covarset().values():
            x_ii = 0
            jj = 0
            for name in covar.var_dp:
                flag = False

                for ii, var in self.opt_prob.get_varset().iteritems():
                    if name == var.name:
                        x_ii += covar.slope_[jj] * self._x[ii]
                        flag = True
                        break

                if not flag:
                    raise ValueError('{} is not in the variable set!'
                                     .format(name))

                jj += 1

            x_ii += covar.intercept_

            xc.append(x_ii)

        self._xc = xc

    def _update_input(self):
        """Update the input file.

        Key words in the template of the input file which will be
        substituted should be put between '<' and '>'.
        """
        with open(self.input_file, 'w') as fp:
            ptns = set()
            for line in self._input_template_lines:
                loop = 0
                while loop < 100:
                    loop += 1

                    # Comment line starting with '!'
                    if re.match(r'^\s*!', line):
                        break

                    # Cannot find '<' or '>'
                    if line.find('<') < 0 or line.find('>') < 0:
                        break

                    # If '<' is on the right of '>'
                    if line.find('<') >= line.find('>'):
                        break

                    # In line comment
                    if line.find('<') > line.find('!') >= 0:
                        break

                    ptn = line[line.find('<')+1:line.find('>')]
                    find_ptn = 0
                    for ii, var in self.opt_prob.get_varset().iteritems():
                        if find_ptn == 1:
                            break
                        if ptn == var.name:
                            line = line.replace('<' + ptn + '>',
                                                str(self._x[ii]), 1)
                            find_ptn = 1

                    for ii, covar in self.opt_prob.get_covarset().iteritems():
                        if find_ptn == 1:
                            break
                        if ptn == covar.name:
                            line = line.replace('<' + ptn + '>', str(self._xc[ii]), 1)
                            find_ptn = 1

                    for ii, staticvar in self.opt_prob.get_staticvarset().iteritems():
                        if find_ptn == 1:
                            break
                        if ptn == staticvar.name:
                            line = line.replace('<' + ptn + '>', str(staticvar.value), 1)
                            find_ptn = 1

                    if find_ptn == 0:
                        raise ValueError(
                            '{} is not any kind of variables!\n'.format(ptn))
                    else:
                        ptns.add(ptn)

                fp.write(line)

        # Looking for unused variables
        for var in self.opt_prob.get_varset().values():
            if var.name not in ptns:
                raise ValueError("Variable '{}' is not found in the input file!"
                                 .format(var.name))

        for covar in self.opt_prob.get_covarset().values():
            if covar.name not in ptns:
                raise ValueError("Co-variable '{}' is not found in the input file!"
                                 .format(covar.name))

        for staticvar in self.opt_prob.get_staticvarset().values():
            if staticvar.name not in ptns:
                raise ValueError("Static-variable '{}' is not found in the input file!"
                                 .format(staticvar.name))

    def _remove_output_files(self):
        """Remove files generated in each simulation"""
        for value in self.fit_points.__dict__.values():
            if isinstance(value, PhaseSpace):
                try:
                    os.remove(value.particle_file)
                except OSError:
                    pass

        for value in self.sections.__dict__.values():
            if isinstance(value, BeamEvolution):
                if value.particle_type == 'astra':
                    for suffix in ['.Xemit.001', '.Yemit.001', '.Zemit.001', '.TRemit.001']:
                        try:
                            os.remove(value.root_name + suffix)
                        except OSError:
                            pass
                elif value.particle_type == 'impact':
                    for suffix in ['.18', '.24', '.25', '.26']:
                        print value.root_name
                        try:
                            os.remove(value.root_name + suffix)
                        except OSError:
                            pass

    def _run_simulation(self):
        """Run the simulation."""
        self._remove_output_files()

        if self.complete_shell is True:
            subprocess.call(self.run_code, shell=True)
        else:
            subprocess.call("timeout {}s {} {}>/dev/null".format(
                self.time_out, self.run_code, self.input_file), shell=True)

    def _write_history(self, dt):
        """Save the optimization history to a log file."""
        if self._n_iter == 1:
            with open(self.log_file, 'wb') as fp:
                fp.write(self.opt_prob.name
                         + ', optimizer = {}'.format(self._optimizer_name)
                         + ', nrestart = {}'.format(self._n_restart)
                         + '\n\n')

                fp.write('{:>6}, {:>12}'.format('niter', 'runtime (s)'))

                for obj in self.opt_prob.get_objset().values():
                    fp.write(", {:>16}".format(obj.name))

                for con in self.opt_prob.get_conset().values():
                    fp.write(", {:>16}".format(con.name))

                for var in self.opt_prob.get_varset().values():
                    fp.write(", {:>16}".format(var.name))

                for covar in self.opt_prob.get_covarset().values():
                    fp.write(", {:>16}".format(covar.name))

                for staticvar in self.opt_prob.get_staticvarset().values():
                    fp.write(", {:>16}".format(staticvar.name))

                fp.write('\n')

        with open(self.log_file, 'a') as fp:
            fp.write('{:6d}, {:12.4f}'.format(self._n_iter, dt))

            for ele in self._f:
                fp.write(', {:16.6e}'.format(ele))

            for ele in self._g:
                fp.write(', {:16.6e}'.format(ele))

            for ele in self._x:
                fp.write(', {:16.6e}'.format(ele))

            for ele in self._xc:
                fp.write(', {:16.6e}'.format(ele))

            if self.opt_prob.get_staticvarset().keys():
                staticvar_values = [staticvar.value for staticvar in
                                    self.opt_prob.get_staticvarset().values()]
                for ele in staticvar_values:
                    fp.write(', {:16.6e}'.format(ele))

            fp.write('\n')

    def _solution_text(self):
        """Write the solution to file"""
        ftext = "\nSolution with optimizer pyOpt.{} on \n{} {}\n".format(
            self._optimizer_name.upper(), self.run_code, self.input_file)
        ftext += "Total time - {}".format(self._time_consumption)
        ftext += self.opt_prob.__str__()
        ftext += self.fit_points.__str__()
        ftext += self.sections.__str__()

        return ftext


def quad_k2g(k, p):
    """Convert the K value of a quadrupole to gradient.

    Parameters
    ----------
    k: float
        Quadrupole strength (1/m^2)
    p: float
        Normalized momentum

    Returns
    -------
    Quadrupole gradient (T/m)
    """
    me = 9.10938291e-31
    v_light = 299792458
    qe = 1.60217657e-19

    return -1.0*k*p*me*v_light/qe
