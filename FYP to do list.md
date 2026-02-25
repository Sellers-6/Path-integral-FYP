\- Implement a multithreaded option to the C++ code, consider adding "#pragma omp simd" to forcibly vectorise certain loops


\- Read section in Creutz Freedman that details the numerics of the DWP - Mackenzie only included tunneling up to an unknown constant which was non-trivial to calculate.


\- David grabovsky Instanton paper


\- WKB understanding - Sakurai, this resource is somewhat limited in its applications to the DWP. We can still use WKB to display the necessity of Numerical estimations in certain systems.
" the Wentzel–Kramers–Brillouin (WKB) method [1], which can be derived as
a semi-classical approximation to the Feynman path integral [20–22], is particularly effec-
tive at producing accurate estimates for the eigenvalues of the quartic oscillator and other
potentials [18]. " - Path integral Monte Carlo method for the quantum anharmonic oscillator


\- Use Dirichlet BCs to enforce these conditions from Mackenzie "As a final note, we have calculated the PI with q = q′ = a; a good exercise is to do the analogous calculation for q = −a, q′ = a". Both of these can be simulated and have their energy splitting calculated.


\- Properly implement and analyse the anharmonic potentials, find their exact solutions using the schroedinger picture. (Use Mittal as a reference)


\- Learn about "Spline fitting" and possibly use this to connect "knots" in the DWP or anharmonic oscillator plots.


\- Think of arguements as to why we are using the path integral rather than the operator/schrodinger methods. One is that we cannot perturb the DWP in schrodinger???


\- Read about how the Lattice size and spacing effects error (it is very small)


\- Write up sections of the report that have been understood


Extensions:


\- Research how spin can be incorporated into instanton/QHO, understand Grassmann variables and incorporate this into the existing C++ code.


\- Observe what changed when incorporating spin and attempt to understand this