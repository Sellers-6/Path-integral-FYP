\- Switch MCMC code over to Rust from C++


\- Read up on theory. This includes:
\- Read all of David grabovsky Instanton lecture notes
\- WKB understanding - Sakurai, this resource is somewhat limited in its applications to the DWP. We can still use WKB to display the necessity of Numerical estimations in certain systems. " the Wentzel–Kramers–Brillouin (WKB) method [1], which can be derived as a semi-classical approximation to the Feynman path integral [20–22], is particularly effective at producing accurate estimates for the eigenvalues of the quartic oscillator and other potentials [18]. " - Path integral Monte Carlo method for the quantum anharmonic oscillator


\- Properly implement and analyse the anharmonic potentials, find their exact solutions using the schroedinger picture. (Use Mittal as a reference)


\- Think of arguements as to why we are using the path integral rather than the operator/schrodinger methods. One is that we cannot perturb the DWP in schrodinger???


\- Read about how the Lattice size and spacing effects error (it is very small)


\- Write up sections of the report that have been understood


Extensions:


\- Research how spin can be incorporated into instanton/QHO, understand Grassmann variables and incorporate this into the existing C++ code.


\- Observe what changed when incorporating spin and attempt to understand this