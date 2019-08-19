"""!
@file Transport.py
@package Coeus

@defgroup Transport Transport

@brief TBD

@author James Bevins

@date 18Aug19
"""

import logging

import numpy as np

module_logger = logging.getLogger('Coeus.MCNPUtilities')
#-----------------------------------------------------------------------------#
class Transport():
    """!
    @ingroup Transport
    The class creates a Transport object to store the transport input file,
    parse the associated sampled variables, and build the transport output.
    """

    ##
    def __init__(self, transInputPath, transCode='mcnp6'):
        """!
        Constructor to build the Transport class.

        @param transInputPath: \e string \n
            The path to the transport input file. \n
        """

        ## @var transInputPath: \e string
        # A path for the transport input file. \n
        self.transPath = transInputPath
        ## @var code: \e string
        # The transport code used. \n
        self.code = transCode

        ## @var sampVars: \e dictionary
        # A dictionary containing the sampled variables. \n
        self.sampVars = {}
        ## @var corrVars: \e dictionary
        # A dictionary containing the correllated variables that are not
        # directly sampled. \n
        self.corrVars = {}
        ## @var transInput: \e string
        # The transport input. \n
        self.transInput = ""
		
		# Populate class variables
        self.importfile()

    def __repr__(self):
        """!
        Transport print function.

        @param self: <em> object pointer </em> \n
            The Transport pointer. \n
        """
        return "Transport({}, {}, {}, {})".format(self.sampVars,
                                                  self.corrVars,
                                                  self.code,
                                                  self.transInput)

    def __str__(self):
        """!
        Human readable Transport print function.

        @param self: <em> object pointer </em>\n
            The Transport pointer. \n
        """

        header = ["Transport:"]
        header += ["Sampled Variables: {}".format(self.sampVars)]
        header += ["Correlated: {}".format(self.corrVars)]
        header += ["Transport Code: {}".format(self.code)]
        header += ["Transport Input: {}".format(self.transInput)]
        return "\n".join(header)+"\n"

    def importfile(self):
        """
        This section imports the first part of the transport input file. It
        breaks it into dictionaries containing sampled variables, correlated
        variables, and the main transport input.

        The key assumption to the import are:
            1) The sampled variables are located at the top of the file and
            are proceded by a '#'
            2) The correlated variables are located next and are proceded by
            a '@'
            3) The variables are separated from the main input by a blank line.

        @param self: <em> object pointer </em>\n
            The Transport pointer. \n
        """
        # Open file
        f = open(self.transPath)

        # Parse input into sampled dict and correlated dict
        for line in f:
            if line[0] == '#':
                line = line.rstrip('\n').replace('# <', '').replace(' =', '=')
                variable, value = line.split('=')
                self.sampVars[variable] = value
            elif line[0] == '@':
                line = line.rstrip('\n').replace('@ <', '').replace(' =', '=')
                variable, value = line.split('=')
                self.corrVars[variable] = value
            elif line in ['\n', '\r\n']:
                break

        # Store main body of input
        self.transInput = f.read()

## Read the generated MCNP output and return the tally results
#  @param path String The path, including filename, to the MCNP output file to be read
#  @param tnum String The number of the tally to be read
#  @return tally array Array of tally results  
def Read_Tally_Output(path, tnum):
    
    assert isinstance(path, str)==True, 'Path must be of type str.'
    assert isinstance(tnum, str)==True, 'tnum must be of type str.'
    
    # Initialize the tally
    tally=[]
    t=False
    
    # Create and open input file 
    try:
        with open(path, "r") as f:
            
            # Read the output file line by line
            for line in f:
            
                # Find key word for start of flux array
                split_list=line.strip().split()
                if len(split_list)>=3:
                    if split_list[0].strip()=="1tally" and split_list[1].strip()==tnum.strip() and split_list[2].strip()=="nps":
                        t=True
                        # Advance 12 lines
                        for i in range(0,11):
                            split_list=f.next().split()
                
                if t==True:
                    # Exit at end of Tally
                    if split_list[0].strip()== "total":
                        t=False
                    else:
                        # Store Tally
                        tally.append([float(split_list[0].strip()),float(split_list[1].strip())])

        # Close the file
        f.close()
    
    except IOError as e:
        module_logger.error("I/O error({0}): {1}".format(e.errno, e.strerror))    
        module_logger.error("File not found was: {0}".format(path))

    # Test that the file closed
    assert f.closed==True, "File ({}) did not close properly.".format(path)
    
    return np.asarray(tally)   

## Read the generated MCNP output and return the tally results
    # @param path String The path, including filename, to the MCNP output file to be read
    # @param tnum String The number of the tally to be read.  Returns the entire binned tally.
    # @param rnum String The number of the tally to be read for the total reactions only. 
    # @return tally array Array of tally results for the tally specified by tnum [Ebins, tally, uncertainty]
    # @return rxs array Total number of reactions for the tally specified by rx_num [tally, uncertainty]
    # @return weight float The total weight of the system
def Read_MCNP_Output(path, tnum, rnum):

    assert isinstance(path, str)==True, 'Path must be of type str.'
    assert isinstance(tnum, str)==True, 'tnum must be of type str.'
    assert isinstance(rnum, str)==True, 'rnum must be of type str.'
    
    # Initialize the tally
    tally=[]
    t=False
    r=False
    w=False
    
    # Create and open input file 
    try:
        with open(path, "r") as f:
            
            # Read the output file line by line
            for line in f:
            
                # Find key word for start of flux array
                split_list=line.strip().split()
                if len(split_list)>=3:
                    if split_list[0].strip()=="1tally" and split_list[1].strip()==tnum.strip() and split_list[2].strip()=="nps":
                        t=True
                        # Advance 11 lines
                        for i in range(0,11):
                            split_list=f.next().split()
                    if split_list[0].strip()=="1tally" and split_list[1].strip()==rnum.strip() and split_list[2].strip()=="nps":
                        r=True
                        # Advance 12 lines
                        for i in range(0,12):
                            split_list=f.next().split()
                    if split_list[0].strip()=="cell" and split_list[1].strip()=="mat" and split_list[2].strip()=="density":
                        w=True
                        # Advance 2 lines
                        for i in range(0,2):
                            split_list=f.next().split()
                
                # Store data if keyword located
                if t==True:
                    # Exit at end of Tally
                    if split_list[0].strip()== "total":
                        t=False
                    else:
                        # Store Tally
                        tally.append([float(split_list[0].strip()),float(split_list[1].strip())])
                        
                if r==True:
                    # Store the total and exit at end of Tally
                    if split_list[0].strip()== "total":
                        rxs=[float(split_list[1].strip()),float(split_list[2].strip())]
                        r=False
                        
                if w==True:
                    # Store the total and exit at end of Tally
                    if len(split_list) > 0:
                        if split_list[0].strip()== "total":
                            weight=float(split_list[2].strip())
                            w=False

        # Close the file
        f.close()
   	 
   	    # Test that the file closed
        assert f.closed==True, "File ({}) did not close properly.".format(path)    

    except IOError as e:
        module_logger.error("I/O error({0}): {1}".format(e.errno, e.strerror))
        module_logger.error("File not found was: {0}".format(path))

    return np.asarray(tally), np.asarray(rxs), weight
