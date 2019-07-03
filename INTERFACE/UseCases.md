###This document outlines possible use cases.  It is not intended to be exaustive, but instead to define some specific examples and how they would interact with the code and current capabilities. The tecnhical nuclear forensics describes the base case that the code was developed for.

##Technical Nuclear Forensics (TNF):
###Overview
For TNF, the goal is to generate target spectra to replicate those from a nuclear weapons to get the correct fission and activation products (debris) for post-detonation attribution.  This case envisioned the National Ignition Facility (NIF) as the starting source, and most of the constraints are derived from fielding at this facility.  The search space has been simplified (2D, no optimization of the envelop, 1D search for foil location) to match this problem.  
###Geometry
The Coeus inputs allow for the geometry to be fully specified in 2D for this case. The user has control over all of the key components for the outer envelop, structural material thickness, foil holder, placement in the target chamber, and snout mount.
###Constraints
The contraints considered are weight and number of fissions (reactions).  These are set through the user inputs.  The Read_MCNP_Output function gets the weight and reactions from a set tally structure.  The weight is implemented as a hard constraint, and the # of fissions is incorporated as a soft contraint in the objective function calculation.  
###Objective Function
The objective function is calculated in the RelativeLeastSquares function.  It is a weighted sum of the squares of the differences between the two spectra.  

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
###Challenges
1) Implement a flexible user geometry input that allows for parts of the problem to be fixed and parts to be variable. 

##Radiation Shields and Converters
###Overview
There are many applications covered under this scope that might be useful. Here we describe two of the more interesting/challenging cases.  These numbers are carried forth throughout this section for clarity.
1) Space shielding design: For satelites, space reactors, and future space travel, light weight, efficient shields are required.  This is one area in nuclear engineering where more advanced optimization techniques have been applied in the past, but many of these problems are still solved by a parametric study.
2) Radiation converter design: It is sometimes useful to convert one type of radiation into another.  
###Geometry
The geometries for these applications tend to be fairly straighforward in the traditional layered concept.  However, matrices of material have been introduced and show some promise for advancing each.  The matrices can be handled explicitly by modeling 3D repeating or stochastic lattices or implicity by modleing homogenous material with a varying composition fractions.  Either is available under Coeus, but which is more appropriate depends on the actual material.  
###Constraints
Constraints for both tend to be mass and thickness.  There can be cost factors as well.  Given the simple geometries, this weight is easy to check before running MCNP, which could enhance efficiency by tossing constrained designs and generating a new design before run time.
###Objective Function
The objective function for 
1) Either a dose threshold or a dose reduction percentage.  This will require tally a flux and implementing a flux to dose conversion tally.  This is all on the user input side.
2) This typically is either a target conversion rate (or percentage) or a generated particle type rate. These are both straighforward tallies and can be implemented within the imput from the user.
###Challenges
1) Ensure that the input reader can handle MCNP universes.
2) May require the use of alternative codes. 

##Medical Isotope Production
###Overview
Medical isotope production is a fairly big business, and it appears that the target designs are optimized by experience, back of the envelope calcuulations, parametric studies.  The basic idea here is that you want to generate as much of isotope X as posssible while avoiding other isotopes that are difficult to seperate from the irradiated sample.  Based on the relative cross-section of the two or more reactions, there will be an optimal beam energy that is obtained by a) tuning the incident neutron beam  or b) degraders.  Since beam is expensive and limited, often a stackup of irradiations is performed at once to generate multiple medical isotopes in sequence as the beam downshifts in energy through the irradiated sample and degraders.  These beams can be either charged particle or neutrons.  
###Geometry
Typically a stackup of layered materials.  It maybe possible that a 2D matrix would also be desirable, primarily for neutron beams.  
###Constraints
Typically there are not difficult constraints here outside of physics.  Space an mass can be a limit depending on the facility.  This may be more of a consideration for neutron beams than charged particles.  You could incorporate maximum bad/good isotope ratios here if the goal is to produce the maximum number of atoms up to some contamination percentage. 
###Objective Function
Objective functions can be the sum of the number of atoms created for each isotope (with or without preference weighting) ot the ratio of bad to good isotopes. It is possible that both of these will be the objective function of interest and need to be combined.  Each calculation is possible to be done inside of a tally, but the combination and weighting needs to be handled elsewhere in post-processing.
###Challenges
1) The ability to run other codes. However, MCNP is probably capable of most of this work and would be the most liekly to be used.
2) The ability to have user specified objective functions made from a combination of a multiple tallies. 

##Fusion Blanket Design
###Overview
Generating sustainable fusion is only half of the problem.  Once that happens, to make it economically viable, the energy must be harvested and new tritium must be bred.  This is typically done in a blanket that surrounds the core.  
###Geometry
These typically are in some layered sphereical or eliptical geometry.  There is typically some level of symmetrry, and the problem can be simplified to 1D or perhaps 2D, to accelerate the calculations.  
###Constraints
There is typically a minimum tritium breeding ratio and energy capture rate.  There could also be a maximum tritium breeding ratio. There may also be some wieght/size constraint, but these tend to be rather large.  
###Objective Function
While there may be a lower bounds for the minimum tritium breeding ratio and energy capture rate, these are typically what you want to maximize.  This can be treated as a multi-objective function, or they could be combined into a single objective (either weighted or unweighted).  
###Challenges
1) Implement materials that can varying weight/atomic fractions that follow some requirements/relationships dictated by materials science. 
2) Multi-objective optimization

##Inverse Problems
###Overview
There are many applications covered under this scope that might be useful.  The general concept is that you have a measurement, and you want to know what the source causing that measurement was.  This could be complicated by an intervening environment that is either known or variable.  
###Geometry
Here the geometry can be anything.  It can also be a fixed value or have variable components.  The source can itself be a variable (energy, distribution, size, location, etc). 
###Constraints
The constraints will vary widely depending on the specific application, but some to consider are masses of materials, ordering of materials, types of sources, and source intensities, among others. 
###Objective Function
This can be formulated many ways but in essence is a comparison to data measurement(s)/observation(s).  This will likely require post-processing of multiple tallies.  It is possible that multi-ojective optimization would be an improvement, but it is believed that these can be solved with single objective optimization.
###Challenges
1) It is likely that a different code will be desired for this application, so that interface needs to be easy to move from MCNP to something else.
2) Ensure variable sources are addressed by the inputs
3) Multi-objective optimization
4) The ability to have user specified objective functions made from a combination of a multiple tallies. 

##Reactor Design
###Overview
These range from assembly level to full core calculations.  
###Geometry
###Constraints
###Objective Function
###Challenges
1) The big question here is if the cost of the calculation can be brought in-line with the number of function evaluations required to optimize. This is not just a question of radiation transport efficiency however.  The choice of the number of variables, variable domain, weight windows, and code will go a long ways to determining the run time.  This is outside of Gnowee/Coues development, but something to consider.
2) It is likely that a different code will be desired for this application, so that interface needs to be easy to move from MCNP to something else.
