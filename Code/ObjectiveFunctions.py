
"""!
@file ObjectiveFunctions.py
@package Coeus

@defgroup ObjectiveFunctions ObjectiveFunctions

@brief Defines a class to perform objective function calculations.

@author James Bevins, Youdong Zhang

@date 22April
"""

_FUNC_DICT = {"relative_least_squares": relative_least_squares,
              "least_squares": least_squares, "u_opt": u_opt}

#-----------------------------------------------------------------------------#    
class ObjectiveFunction:
    """!
    @ingroup ObjectiveFunctions
    The class creates a ObjectiveFunction object that can be used in
    optimization algorithms.
    """

    def __init__(self, method=None, tallyNum=None, objective=None):
        """!
        Constructor to build the ObjectiveFunction class.

        @param self: <em> object pointer </em> \n
            The object pointer. \n
        @param method: \e string \n
            The name of the objective function to evaluate. \n
        @param tallyNum: \e string \n
            An associated MCNP tally number that is to be used to provide the
            input for the objective function calculation. \n
        @param objective: <em> integer, float, or numpy array </em> \n
            The desired objective associated with the optimization.  The
            chosen value and type must be compatible with the optiization
            function chosen. \n
        """

        ## @var funcName <em> function handle </em> The function handle for
        # the objective function to be used for the optimization.  The
        # function must be specified as a method of the class.
        self.funcName = method 
        ## @var funcTally \e integer The MCNP tally number to be used to
        #provide the input for the objective function calculation.
        self.funcTally = tallyNum 
        ## @var objective  <em> integer, float, or numpy array </em> The
        # desired outcome of the optimization.
        self.objective = objective 

    def __repr__(self):
        """!
        ObjectiveFunction class param print function.

        @param self: \e pointer \n
            The ObjectiveFunction pointer. \n
        """
        return "ObjectiveFunction({}, {})".format(self.funcName,
                                                  self.funcTally)

    def __str__(self):
        """!
        Human readable ObjectiveFunction print function.

        @param self: \e pointer \n
            The ObjectiveFunction pointer. \n
        """

        header = ["\ObjectiveFunction:"]
        header += ["Objective Function: {}".format(funcName._Name_)]
        header += ["Tally Number: {}".format(funcTally)]
        return "\n".join(header)+"\n"
 
    def Uopt(self, c):
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

        assert len(c)==len(self.objective), "The length of the candidate and \ 
                                objective  must be equal in Uopt."  

        return np.sum(abs(self.objective-c))

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

        assert len(c)==len(self.objective), "The length of the candidate and \ 
                                objective  must be equal in least_squares."  

        return np.sum((self.objective-c)**2)

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

        assert len(c)==len(self.objective), "The length of the candidate and \ 
                                objective  must be equal in \
                                relative_least_squares." 

        # For bins with no tally results, project the fitness using simple
        # linear extrapolation
        if project == True:
            for i in range(len(c)):
                if c[i]==0.0:
                    extrapIndex1 = i + 1
                    extrapIndex2 = i + 2
                    while c[extrapIndex1] == 0.0 or c[extrapIndex2] == 0.0:
                        extrapIndex1, extrapIndex2 += 1
                    c[i] = c[extrapIndex1]-(extrapIndex1-i)\
                            *(c[extrapIndex2]-c[extrapIndex1]
                              /(extrapIndex2-extrapIndex1))
        return np.sum((self.objective-c)**2/self.objective)