**What the code does:**
It calculates the electrostatic potential in a space bounded by two charged plates. The space between the plates is half filled with liquid 1 and the other half is filled with liquid 2 and the plates carry constant surface charge densities.

**How does it do so:**
It minimizes the relevant free energy by using a gradient descent technique. The solutions are effectively the solution of a nonlinear differential equation (Poisson-Boltzmann equation) where we expand the nonlinear term using a series expansion and the parameter 'n' in the code controls the degree of nonlinearity we allow for. Thus, n=1 corresponds to a linear problem (also known as Debye-Huckel problem) and n->infinity corresponds to the full nonlinear (Poisson-Boltzmann problem). However, it turns out that n=4 is already sufficient as the differences between the solutions for n=4 and n=5 are negligible for most practical purposes. For further details please have a look at our published paper: [J. Chem. Phys. **145**, 064707 (2016)](https://doi.org/10.1063/1.4960623) | [arXiv:1606.03908](https://doi.org/10.48550/arXiv.1606.03908)

**Sample parameters and input files:**
After compiling, the code can be run with the parameters in the runscript.txt file. Expected output files are also provided.

At the left and right boundaries of the system, i.e., at the two surfaces, constant surface charge density boundary condition is used. However, for running the code, one also needs two files file_phi_low.txt and file_phi_up.txt 
providing values of the potentials at the lower and upper boundaries, respectively. They are provided. The source codes to obtain them are also provided in the subfolders.

For n=1, the problem is linear and exact analytical solutions are available. Those can be computed for each setup using the Mathematica scripts named 'Exact.nb' (or using scripts under the subfolders named 'Exact_potentials') to 
check the accuracy. Since this is a linear problem, one can also solve it numerically using an alternate technique. The energy minimization results in a system of linear equations which can be solved using Gauss-Seidel method. 
The corresponding code and sample output files are provided for comparison inside the sobfolder named 'Debye_Huckel'.
