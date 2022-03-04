# BARW
Includes scripts for simulations of Branching and annihilating random walks (BARWs) that are
guided by an external potential and can additionally have local self-interactions. 
Theory and technical details are described in the publication:

Ucar et al, Nat Commun. 12, 6830 (2021)
https://www.nature.com/articles/s41467-021-27135-5

- The script branching_rules.py contains the rules underlying the BARW framework,
as well as additional modules for the implementation of the external guidance and 
local self-interaction effects.

- The script branching_simulation.py defines a class for branched networks, and 
includes a module to generate simulation loops.

- The script test_simulation.py provides a simple test case with some selected
network properties such as its branching rate, presence of external guidance etc. 
These parameters can be changed to explore different scenarios.
