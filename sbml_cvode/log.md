# Version #

version 1.1

# Todo #

*	Units to .xml and Param.cpp file

# Development history #

## version 1 ##

1.	Reactions
	-	Reactions need to have kinetic law defined.
2.	Events
	-	Event trigger
	-	Event assignment
		-	(Delay is not yet fully supported)
3.	Initial Assignments
4.	Rules:
	-	Assginment rules are supported
5.	Parameters:
	-	only global parameters

### version 1.1 ###

*	(Static) class parameters are no longer const and can be altered after instantiation.
*	External parameter files to set all variable values.
*	Simulation time and tolerance from parameter file
*	Evaluate initial assignment after reading parameter file.
*	Normal and Lognormal parameter sampling; two-stage parameter sampling.
*	Serialization (including class parameters) tested.

