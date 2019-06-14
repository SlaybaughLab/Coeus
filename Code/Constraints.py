"""!
@file Constraints.py
@package Coeus

@defgroup Constraints Constraints

@brief Defines a class to perform constraint calculations.

@author James Bevins, Youdong Zhang

@date 14Jun19
"""

import logging

from math import ceil

module_logger = logging.getLogger('Coeus.Constraints')

#-----------------------------------------------------------------------------#
class Constraints(object):
    """!
    @ingroup Constraints
    The class creates a Constraints object that can be used in
    optimization algorithms.
    """

    def __init__(self, method=None, constraint=None, tallyNum=None,
                 penalty=1E15):
        """!
        Constructor to build the ObjectiveFunction class.

        @param self: <em> object pointer </em> \n
            The object pointer. \n
        @param _FUNCT_DICT: \e dictionary \n
            A mapping from string function names to function handles. \n
        @param method: \e string \n
            The name of the constraint function to evaluate. \n
        @param constraint: \e float \n
            The constraint to be compared against. \n
        @param tallyNum: \e integer \n
            The tally number associated with the constraint. \n
        @param penalty: \e float \n
            The penalty to be applied if a constraint is violated.  1E15
            is recommended. \n
        """

        ## @var _FUNC_DICT <em> dictionary of function handles </em> Stores
        # the mapping between the string names and function handles for
        # the constraint function evaluations in the class.  This must be
        # updated by the user if a function is added to the class.
        self._FUNC_DICT = {"less_or_equal":
                           self.less_or_equal,
                           "greater_than": self.greater_than}
        ## @var func <em> function handle </em> The function handle for
        # the constraint function to be used for the optimization.  The
        # function must be specified as a method of the class.
        if method is not None:
            self.func = self.set_constraint_func(method)
        else:
            self.func = method
        ## @var constraint \e float The constraint to be enforced.
        self.constraint = constraint
        ## @var tallyNum \e float The tally associated with the constraint.
        self.tallyNum = tallyNum
        ## @var penalty \e float The penalty to be applied if the constraint
        # is violated
        self.penalty = penalty

        module_logger.info('User defined inputs: {}'.format(print(self)))

    def __repr__(self):
        """!
        Constraint class param print function.

        @param self: \e pointer \n
            The Constraint pointer. \n
        """
        return "Constraint({}, {}, {}, {})".format(self.func.__name__,
                                              self.constraint,
                                              self.tallyNum,
                                              self.penalty)

    def __str__(self):
        """!
        Human readable Constraint print function.

        @param self: \e pointer \n
            The Constraint pointer. \n
        """

        header = ["Constraint:"]
        header += ["Constraint Function: {}".format(self.func.__name__)]
        header += ["Constraint: {}".format(self.constraint)]
        header += ["Constraint Tally: {}".format(self.tallyNum)]
        header += ["Penalty: {}".format(self.penalty)]
        return "\n".join(header)+"\n"

    def set_constraint_func(self, funcName):
        """!
        Converts an input string name for a function to a function handle.

        @param self: \e pointer \n
            The Constraint pointer. \n
        @param funcName \e string \n
             A string identifying the constraint function to be used. \n
        """
        self.func = self._FUNC_DICT[funcName]
        assert hasattr(self.func, '__call__'), 'Invalid function handle'

    def get_penalty(self, violation):
        """!
        Calculate the constraint violation penalty, if any.

        @param self: \e pointer \n
            The Constraint pointer. \n
        @param violation \e float \n
             The magnitude of the constraint violation used for scaling the
             penalty. \n

        @return \e float: The scaled penalty. \n
        """

        return self.penalty*ceil(violation)**2

#-----------------------------------------------------------------------------#
# The following sections are user modifiable to all for the use of new
# objective functions that have not yet been implemented.  The same format must
# be followed to work with the standard Coeus call. If a function is added.
# it must also be added to the _FUNC_DICT attribute of the class.
#-----------------------------------------------------------------------------#

    def less_or_equal(self, candidate):
        """!
        Compares a previously calculated value to a user specifed maximum.

        @param self: \e pointer \n
            The Constraint pointer. \n
        @param candidate \e float \n
             The calculated value corresponding to a candidate design.\n

        @return \e float: The penalty associated with the candidate design. \n
        """

        if candidate <= self.constraint:
            return 0

        return self.get_penalty(candidate-self.constraint)

    def greater_than(self, candidate):
        """!
        Compares the calculated value to the minimum specified by the user.

        @param self: \e pointer \n
            The Constraint pointer. \n
        @param candidate \e float \n
             The calculated value corresponding to a candidate design.\n

        @return \e float: The penalty associated with the candidate design. \n
        """

        if candidate > self.constraint:
            return 0

        return self.get_penalty(self.constraint - candidate)
