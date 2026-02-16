\- Implement a multithreaded option to the C++ code, consider adding "#pragma omp simd" to forcibly vectorise certain loops


\- Remove the correlation function vector from C++. This can be evaluated from the position data so can be done locally in R. It currently makes file size ~ 2x bigger


\- WKB understanding - Shculmann and Sakurai, use this to evaluate accuracy of simulated values


\- Read about how the Lattice size and spacing effects error (it is very small)


\- Write up sections of the report that have been understood


Extensions:


\- Research how spin can be incorporated into instanton/QHO, understand Grassmann variables and incorporate this into the existing C++ code.


\- Observe what changed when incorporating spin and attempt to understand this