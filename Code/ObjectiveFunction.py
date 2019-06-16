"""!
@file ObjectiveFunction.py
@package Coeus

@defgroup ObjectiveFunction ObjectiveFunction

@brief Defines a class to perform objective function calculations.

@author James Bevins, Youdong Zhang

@date 14June19
"""

import logging

import numpy as np

module_logger = logging.getLogger('Coeus.ObjectiveFunction')

#-----------------------------------------------------------------------------#
class ObjectiveFunction(object):
    """!
    @ingroup ObjectiveFunction
    The class creates a ObjectiveFunction object that can be used in
    optimization algorithms.
    """

    def __init__(self, method=None, tallyNum=None, objType=None,
                 objForm=None, objective=None):
        """!
        Constructor to build the ObjectiveFunction class.

        @param self: <em> object pointer </em> \n
            The object pointer. \n
        @param _FUNCT_DICT: \e dictionary \n
            A mapping from string function names to function handles. \n
        @param method: \e string \n
            The name of the objective function to evaluate. \n
        @param tallyNum: \e string \n
            An associated MCNP tally number that is to be used to provide the
            input for the objective function calculation. \n
        @param objType: \e string \n
            The type of objective.  Valid entries are "spectrum" or
-           "total". \n
        @param objForm: \e integer \n
            The type of objective.  Valid entries are 0-4. \n
            0 = "mcnp" \n
            1 = "normalized" \n
            2 = "differential" \n
            3 = "normalized_differential" \n
            4 = "lethargy" \n
            5 = "normalized_lethargy" \n
        @param objective: <em> integer, float, or numpy array </em> \n
            The desired objective associated with the optimization.  The
            chosen value and type must be compatible with the optiization
            function chosen. \n
        """

        ## @var _FUNC_DICT <em> dictionary of function handles </em> Stores
        # the mapping between the string names and function handles for
        # the objective function evaluations in the class.  This must be
        # updated by the user if a function is added to the class.
        self._FUNC_DICT = {"relative_least_squares":
                           self.relative_least_squares,
                           "least_squares": self.least_squares,
                           "u_opt": self.u_opt}
        ## @var func <em> function handle </em> The function handle for
        # the objective function to be used for the optimization.  The
        # function must be specified as a method of the class.
        if method is not None:
            self.set_obj_func(method)
        else:
            self.func = method
        ## @var funcTally \e string The MCNP tally number to be used to
        #provide the input for the objective function calculation.
        self.funcTally = tallyNum
        ## @var objType \e string The type of objective function calculation.
        self.objType = objType
        ## @var objForm \e string The form of objective function.  Only
        # specified if the objType is "spectrum".
        self.objForm = objForm
        ## @var objective  <em> integer, float, or numpy array </em> The
        # desired outcome of the optimization.
        self.objective = objective

    def __repr__(self):
        """!
        ObjectiveFunction class param print function.

        @param self: \e pointer \n
            The ObjectiveFunction pointer. \n
        """
        return "ObjectiveFunction({}, {}, {}, {}, {})".format(
                                                         self.func.__name__,
                                                         self.funcTally,
                                                         self.objType,
                                                         self.objForm,
                                                         self.objective)

    def __str__(self):
        """!
        Human readable ObjectiveFunction print function.

        @param self: \e pointer \n
            The ObjectiveFunction pointer. \n
        """

        header = ["ObjectiveFunction:"]
        header += ["Objective Function: {}".format(self.func.__name__)]
        header += ["Tally Number: {}".format(self.funcTally)]
        header += ["Objective Function Type: {}".format(self.objType)]
        header += ["Objective Form: {}".format(self.objForm)]
        header += ["Objective: {}".format(self.objective)]
        return "\n".join(header)+"\n"

    def set_obj_func(self, funcName):
        """!
        Converts an input string name for a function to a function handle.

        @param self: \e pointer \n
            The ObjectiveFunction pointer. \n
        @param funcName \e string \n
             A string identifying the objective function to be used. \n
        """
        self.func = self._FUNC_DICT[funcName]
        assert hasattr(self.func, '__call__'), 'Invalid function handle'

#-----------------------------------------------------------------------------#
# The following sections are user modifiable to all for the use of new
# objective functions that have not yet been implemented.  The same format must
# be followed to work with the standard Coeus call. If a function is added.
# it must also be added to the _FUNC_DICT attribute of the class on line 62.
#-----------------------------------------------------------------------------#

    def u_opt(self, c):
        """!
        Calculated the fitness of a series of values using the U-Optimality
        condition.  From: "Relationships among Several Optimality Criteria
        E" by E.A. Rady.

        @param self: \e pointer \n
            The ObjectiveFunction pointer. \n
        @param c <em> numpy array </em> \n
             The array of results corresponding to a candidate design.  For
             example, this can be an energy spectra of a flux. \n

        @return \e float: The u-optimality criteria based fitness for a
             design. \n
        """

        assert len(c) == len(self.objective), ("The length of the candidate "
                                "and objective  must be equal in u_opt.")

        return np.sum(abs(self.objective[:, 1]-c))

    def least_squares(self, c):
        """!
        Calculated the fitness of a series of values using the least squares
        condition.

        @param self: \e pointer \n
            The ObjectiveFunction pointer. \n
        @param c <em> numpy array </em> \n
             The array of results corresponding to a candidate design.  For
             example, this can be an energy spectra of a flux. \n

        @return \e float: The least squares criteria based fitness for a
             design. \n
        """

        assert len(c) == len(self.objective), ("The length of the candidate "
                              "and objective  must be equal in least_squares.")

        return np.sum((self.objective[:, 1]-c)**2)

    def relative_least_squares(self, c, project=True):
        """!
        Calculated the fitness of a series of values using the relative least
        squares condition. This provides equal weighting to all bins in the
        data set being evaluated, regardless of overall magnitude.

        @param self: \e pointer \n
            The ObjectiveFunction pointer. \n
        @param c <em> numpy array </em> \n
             The array of results corresponding to a candidate design.  For
             example, this can be an energy spectra of a flux. \n
        @param project \e boolean \n
             A flag on wether to project a reasonable guess to bins that have
             zero for values. The projected value is a simple linear
             exprapolation. \n

        @return \e float: The relative_least_squares criteria based fitness for
            a design. \n
        """

        assert len(c) == len(self.objective), ("The length of the candidate "
                      "and objective must be equal in relative_least_squares.")

        # For bins with no tally results, project the fitness using simple
        # linear extrapolation
        if project:
            for i in range(len(c)):
                if c[i] == 0.0:
                    module_logger.warning('User defined tally contains bins '
					                'with zero counts')
                    extrapIndex1 = i + 1
                    extrapIndex2 = i + 2
                    if extrapIndex2 < len(c):
                        while c[extrapIndex1] == 0.0 or c[extrapIndex2] == 0.0:
                            extrapIndex1 += 1
                            extrapIndex2 += 1
                            if extrapIndex2 >= len(c):
                                extrapIndex1 = i - 2
                                extrapIndex2 = i - 1
                                break
                    else:
                        extrapIndex1 = i - 2
                        extrapIndex2 = i - 1
                    c[i] = c[extrapIndex1]-(extrapIndex1-i)\
                            *(c[extrapIndex2]-c[extrapIndex1]
                              /(extrapIndex2-extrapIndex1))
        return np.sum(((self.objective[:, 1]-c)/self.objective[:, 1])**2\
                      *self.objective[:, 1]/sum(self.objective[:, 1]))
