# Machine-Learning-algorithms-in-quantum-biology

Description:
A supervised machine learning algorithm was used to study an outstanding problem in quantum biology pertaining to the formation of delocalised excitonic wavefunctions across chromophores in photosynthesis. Source: Quantum Dynamics and Decoherence in Photosynthetic Systems (Masters thesis by Leonard Ruocco at the University of Sussex,UK)
Sample MATLAB code is included that runs the ML algorithm on a 3x3 Hamiltonian to replicate the results of a 7x7 known Hamiltonian. 
Therefore this represents a supervised ML optimisation technique that uses data produced by solving the Lindblad quantum master equation
for the 7x7 Hamiltonian as supervision for the optimisation routine. The routine samples the parameter space of the 3x3 Hamiltonian, solving
the quantum master equation each time to calculate the reduced densit matrix and test the results against the 7x7 Hamiltonian results. With
each iteration the average difference of the density matrix for the 3x3 matrix and the 7x7 matrix is minimised serving as the optimisation
target parameter. 

An example output comparing the reduced density matrix of the 3x3 Hamiltonian to the 7x7 Hamiltonian is included in the repository.

A plot showing multiple runs of the ML algorithm as a function of initial energy scaling of the test 3x3 Hamiltonian. This shows the 
existence of local minima in the parameter space.
