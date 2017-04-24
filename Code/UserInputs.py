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

OBJECTIVE FUCNTION PARAMETERS

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
    5 = "normalized lethargy"
    
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

@date 23April
"""

from ObjectiveFunction import ObjectiveFunction

#------------------------------------------------------------------------------#
class UserInputs(object):
    """!
    @ingroup UserInputs
    The class creates a UserInputs object to store the user input file
    locations, read the user inputs, and set the appropriate classes required
    to run Coeus.
    """

    ##
    def __init__(self, coeusInputPath=None, mcnpInputPath=None):
        """!
        Constructor to build the UserInputs class. If paths is specified, the
        object attributes are populated.

        @param self: <em> object pointer </em> \n
            The objeUserInputsct pointer. \n
        @param coeusInputPath: \e string \n
            The path to the coues input file. \n
        @param mcnpInputPath: \e string \n
            The path to the mcnp input file. \n
        """

        ## @var coeusInput: \e string
        # A path for the Coeus input file.
        self.coeusInput = coeusInputPath
        ## @var mcnpInputPath: \e string
        # A path for the MCNP input file.
        self.mcnpInput = mcnpInputPath

    def __repr__(self):
        """!
        UserInputs print function.

        @param self: <em> object pointer </em> \n
            The UserInputs pointer. \n
        """
        return "UserInputs({}, {})".format(self.coeusInput, self.mcnpInput)

    def __str__(self):
        """!
        Human readable UserInputs print function.

        @param self: <em> object pointer </em>\n
            The UserInputs pointer. \n
        """

        header = ["\UserInputs:"]
        header += ["Coeus Input File Path: {}".format(self.coeusInput)]
        header += ["MCNP Input File Path: {}".format(self.mcnpInput)]
        return "\n".join(header)+"\n"

    def read_coeus_settings(self):
        """!
        Reads the input file and creates the corresponding objects and populates
        their attributes.

        @param self: <em> object pointer </em> \n
            The UserInputs pointer. \n

        @return <em> Objective Function Object </em>: An ObjectiveFunction 
            object initialized with the user input parameters. \n
        """
        
        # Create the relevant objects
        objSet = ObjectiveFunctions()
        # Open file
        try: 
            f = open(self.coeusInput, 'r') 
            
            # Read the file line by line and store the values in the ETA_Params object
            for line in f:
                splitList=line.strip().split(" ")
                for case in Switch(splitList[0].strip().lower()):
                    if case('OBJECTIVE FUNCTION PARAMETERS'.lower()):
                        while True:
                            splitList = f.next().strip().split()
                            if case('function'.lower()):
                                objSet.set_obj_func(splitList[1].strip())
                            if case('tally'.lower()):
                                objSet.funcTally = splitList[1].strip()
                            if case('type'.lower()):
                                objSet.objType = splitList[1].strip()
                            if case('objective'.lower()):
                                num = int(splitList[1].strip())
                                objSet.objForm = int(splitList[2].strip())
                                tmp =[]
                                while len(tmp) < num:
                                    splitList = f.next().strip().split()
                                    for i in range(0,len(splitList),2):
                                        tmp.append(float(splitList[i].strip()),
                                                   float(splitList[i+1].strip()))
                                objSet.spectrum = tmp
                            if case():
                                module_logger.warning("A unkown user input was \
                                            found: {} ".format(splitList[0].strip()))
                    if case():
                        module_logger.warning("A unkown section was specified: \
                        {}".format(splitList[0].strip()))
        
            # Close the file
            self.f.close()
        except IOError as e:
            module_logger.error("I/O error({0}): {1}".format(e.errno, e.strerror))    
       
        # Test that the file closed
        assert self.f.closed==True, "File did not close properly."
        
        return objSet
