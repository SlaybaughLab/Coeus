##Boron Neutron Capture Therapy(BNCT)

###Overview
BNCT is a method to irradiate tumors in the brain.  The goal here is to maximize dose to the tumor and minimize dose to the patient. This problem can be attacked in two ways.  
1) You can assume you know the ideal spectrum based on other calculations. In this case, you are trying to optimize the spectrum out the back end of the ETA. 
2) You can treat it as an integral problem and include the model of the patients head.  In this case, you are trying to maximize the dose to tumor to patient ratio.  This is the more robust approach, and may not be that much more difficult to implement.  
###Geometry
The BNCT geometry is similar to the TNF - a  combination of a conical and cylindrical geometry.  However, there are 3-D components that might be worth adding and optimizing.  To do so, the geometry would need to move to 3-D, and the envelop parameters (at least some of them) would need to be optimizable.  
The question of the diagnostic volume is a tricky one to answer with the current code.  The best approach is to tally the spectrum coming out the back of the ETA, but this does require some code modifications to make work.  It also requires an effective "zeroing out" of the current diagnostic volume.  The diagnostic volume could be generalized such that the patient's head could be modeled, but this would require some rethinking of what user parameters are specified and changes to the code (mainly initialization).  It could be done with the current 2-D search parameters.  
###Constraints
The constraints here are patient irradiation time, dose to tumor, and patient dose.  
####Irradiation Time/Dose to Tumor
The irradiation time is a function of the total dose required on the tumor, so this can be fomulated as a minimum number of reactions, given a specified source strength, in a tally for case #2, which is allowable under the current inputs.  
For case #1, much more knowledge about the problem would be needed to implement this constraint. It would be similar to #2, but you are only tallying the flux, so you have to specify the minimum flux.  That is limited because the reaction rate will vary based on the Energy dependent flux, so you may slightly over or underestimate the time.  
The best way to do this is to use case #2.  The current approach of a "diagnostic volume" could be generalized to allow for a patient's head to be specified.  Then, an ability to have user specified tallies (tied to a given constraint or objective function) needs to be added.  Finally, the general time constraint function needs to be added that takes the user specified tally and src strength to compute the time to meet the minimum dose to tumor.  
####Patient Dose
The current code would require you to go in to the code and hardcode a tally that converts the flux to a dose at the back end.  This is only a crude approximation, and there isn't a good way to give the user a lot of control from the input files.
The best way to do this is to use case #2.  The current approach of a "diagnostic volume" could be generalized to allow for a patient's head to be specified.  Then, an ability to have user specified tallies (tied to a given constraint or objective function) needs to be added.  Finally, the general time constraint function needs to be added that takes the user specified tally and src strength to compute the time to meet the minimum dose to tumor.  
###Objective Function
The objective function is to maximize the dose to tumor to dose to patient ratio.  
The current code can do this by taking the known spectrum that accomplishes this minimization as the objective and using the current spectra comparing objective functions. 
The better approach would be to use case #2.  The current approach of a "diagnostic volume" could be generalized to allow for a patient's head to be specified.  Then, an ability to have user specified tallies needs to be added.  The objective function would then look to maximize the ratio of the two doses that were computed in the tallies.    
