# Finite-Element-Method-Euler-Beam
This program shows FEM in action for given Euler beam, initially beam is split into 20 sections (elements) which is 20 + 1 nodes, on specific nodes there are given moment M, force P, distributed load Q, 3 fixed support nodes and on left side beam is rigid fixed into a wall and on right side there is a spring K. 
Only difference between program versions is that beam0.c uses parameter S to split each of 20 beam elements into S more resulting 20 * S elements and 20 * S  + 1 nodes
which will show how hermite polynomials can approximate resulting M and Q
