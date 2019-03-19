# Implementation corresponding to "Advanced Algebraic Attack on Trivium" and "Algebraic Fault Analysis on Katan"

## Solver for quadratic systems with one solution for cryptoanalytic purposes
Needs M4RI. Described in "Advanced Algebraic Attack on Trivium" (https://link.springer.com/chapter/10.1007/978-3-319-32859-1_23).
## Quadratic equation system for Trivium and Katan(32)
Defined in www.ecrypt.eu.org/stream/papersdir/2006/021.pdf (Trivium) and www.cs.technion.ac.il/~orrd/KATAN/CHES2009.pdf (Katan)
## Solve the systems from Katan and Trivium for reduced rounds with algebraic and fault analysis
The specifics of the analyses are given in "Advanced Algebraic Attack on Trivium" (https://link.springer.com/chapter/10.1007/978-3-319-32859-1_23) and "Algebraic Fault Analysis on Katan" (https://eprint.iacr.org/2014/954.pdf) 



## Setup
Implementation is done in C++ and Sage (CAS based on Python enhancements - see http://www.sagemath.org/) 
1. Setup M4RI and compile quadSolver
2. Setup Sage and PolyBoRi. For Sat-Solver inclusion the implementation uses CryptoMiniSat.
3. For flexibility write a configTrivSolver.sage with the path to eqSolver.sage
4. Insert "needed" (depending on your need see trivSolver.sage, katanSolverFull.sage or katanFaultySolver.sage) files to a folder in your path
5. Make yourself familiar to the parser (see trivSolver.sage, katanSolverFull.sage or katanFaultySolver.sage)
6. Toy with given parameters
7. If you encounter problems write to frank.quedenfeld(at)googlemail.com
