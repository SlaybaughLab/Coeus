"""!
@file UserInputs.py
@package Coeus

@defgroup UserInputs UserInputs

@brief Defines a class to perform objective function calculations.

The keyword parameters discussed here are not case sensitive, and are separated
by any number of spaces.  Any other separator will result in the file not
reading correctly. Only allowed key words will be read, and misspelled key words
will be ignored.  Depending on the input being specified, this may or may not
result in a failed run.  Warnings will be issued for unkown user input lines.

The Coeus user inputs will is constructed as follows:

There is one valid section name.  The section name must be:

OBJECTIVE FUNCTION PARAMETERS

The allowed inputs for the OBJECTIVE FUCNTION PARAMETERS section are:

\e function Expects a function name corresponding to a function in the
    ObjectiveFunctions class.  \n
\e tally Expects an integer number corresponding to an MCNP tally. \n
\e type This specifies the type of tally to be used for the objective function
    evaluation.  Valid inputs are "total" or "spectrum".  The actual tally can
    be have more components, but only the specified porion will be used for the
    objective fucntion evaluation. \n
\e objective If the type specified was "total", then this is a float or integer
    representing the desired objective.  If "spectrum" was specified for the
    type, this is the number of bins for the objective spectrum followed by the
    form of the objective spectrum.  The form responses are:

    0 = "mcnp"
    1 = "normalized"
    2 = "differential"
    3 = "normalized_differential"
    4 = "lethargy"
    5 = "normalized_lethargy"

    Based on the response, the proper calculations will be performed to
    translate the mcnp tally to the correct form for objective function
    calculation. The spectrum is specified on the following line. \n

    The spectrum is enetered in the form of energy amount seperated by a single
    (or multiple spaces. A couple of spectrum examples: \n

    objective 4
    1 0.25 2 0.5 5 0.2 10 0.05

    or

    objective 4
    1 0.25
    2 0.5
    5 0.2
    10 0.05

    NOTE: The objective spectrum specified needs to be the same structure that
    is used used for the tally number specified.

@author James Bevins

@date 16June19
"""

import logging

import numpy as np

from ObjectiveFunction import ObjectiveFunction
from Utilities import Switch

module_logger = logging.getLogger('Coeus.UserInputs')

#-----------------------------------------------------------------------------#
class UserInputs(object):
    """!
    @ingroup UserInputs
    The class creates a UserInputs object to store the user input file
    locations, read the user inputs, and set the appropriate classes required
    to run Coeus.
    """

    ##
    def __init__(self, coeusInputPath=None, transInputPath=None,
                 advantgInputPath=None):
        """!
        Constructor to build the UserInputs class. If paths is specified, the
        object attributes are populated.

        @param self: <em> object pointer </em> \n
            The objeUserInputsct pointer. \n
        @param coeusInputPath: \e string \n
            The path to the coues input file. \n
        @param transInputPath: \e string \n
            The path to the trans input file. \n
        @param advantgInputPath: \e string \n
            The path to the advantg input file. \n
        """

        ## @var coeusInput: \e string
        # A path for the Coeus input file. \n
        self.coeusInput = coeusInputPath
        ## @var transInputPath: \e string
        # A path for the transport input file. \n
        self.transInput = transInputPath
        ## @var advantgInputPath: \e string
        # A path for the transport input file. \n
        self.advantgInput = advantgInputPath

        ## @var code: \e string
        # The transport code used. \n
        self.code = 'mcnp6'

    def __repr__(self):
        """!
        UserInputs print function.

        @param self: <em> object pointer </em> \n
            The UserInputs pointer. \n
        """
        return "UserInputs({}, {}, {}, {})".format(self.coeusInput,
                                                   self.transInput,
                                                   self.advantgInput,
                                                   self.code)

    def __str__(self):
        """!
        Human readable UserInputs print function.

        @param self: <em> object pointer </em>\n
            The UserInputs pointer. \n
        """

        header = ["UserInputs:"]
        header += ["Coeus Input File Path: {}".format(self.coeusInput)]
        header += ["Transport Input File Path: {}".format(self.transInput)]
        header += ["ADVANTG Input File Path: {}".format(self.advantgInput)]
        header += ["Transport Code Used: {}".format(self.code)]
        return "\n".join(header)+"\n"

    def read_coeus_settings(self):
        """!
        Reads the input file and creates the corresponding objects and
        populates their attributes.

        @param self: <em> object pointer </em> \n
            The UserInputs pointer. \n

        @return <em> Objective Function Object </em>: An ObjectiveFunction
            object initialized with the user input parameters. \n
        """

        # Initialize the section headers that are valid:
        sectionHeaders = ['problem definition',
                          'objective function parameters']

        # Create the relevant objects
        objSet = ObjectiveFunction()
        # Open file
        try:
            f = open(self.coeusInput, 'r')

            # Read the file line by line and store the values
            for line in f:
                if line.strip().lower() == 'PROBLEM DEFINITION'.lower():
                    line = f.next().strip().lower()
                    while line not in sectionHeaders:
                        # Stop at end of block
                        try:
                            splitList = line.split()
                        except IndexError:
                            break
                        if not splitList:
                            break

                        for case in Switch(splitList[0].strip().lower()):
                            if case('transport-code'):
                                self.code = splitList[1].strip()
                                break
                            if case():
                                module_logger.warning("Unkown user input "
                                "found: {} ".format(splitList[0].strip()))
                                break

                        # Stop at end of file
                        try:
                            line = f.next().strip().lower()
                        except StopIteration:
                            break

                if line.strip().lower() == \
                   'OBJECTIVE FUNCTION PARAMETERS'.lower():
                    line = f.next().strip().lower()
                    while line not in sectionHeaders:
                        # Stop at end of file
                        try:
                            splitList = line.split()
                        except IndexError:
                            break
                        if not splitList:
                            break

                        for case in Switch(splitList[0].strip().lower()):
                            if case('function'):
                                objSet.set_obj_func(splitList[1].strip())
                                break
                            if case('tally'):
                                objSet.funcTally = splitList[1].strip()
                                break
                            if case('type'):
                                objSet.objType = splitList[1].strip()
                                break
                            if case('objective'):
                                num = int(splitList[1].strip())
                                objSet.objForm = int(splitList[2].strip())
                                tmp = []
                                while len(tmp) < num:
                                    splitList = f.next().strip().split()
                                    for i in range(0, len(splitList), 2):
                                        tmp.append([float(splitList[i].strip()),
                                                float(splitList[i+1].strip())])
                                objSet.objective = np.asarray(tmp)
                                break
                            if case():
                                module_logger.warning("Unkown user input "
                                "found: {} ".format(splitList[0].strip()))
                                break

                        # Stop at end of file
                        try:
                            line = f.next().strip().lower()
                        except StopIteration:
                            break
                else:
                    module_logger.warning("A unkown section was specified: "
                                              "{}".format(line.strip()))

            # Close the file
            f.close()
        except IOError as e:
            module_logger.error("I/O error({0}): {1}".format(e.errno,
                                                             e.strerror))

        # Test that the file closed
        assert f.closed, "File did not close properly."

        module_logger.info('The Objective Function: {}'.format(objSet))

        return objSet
