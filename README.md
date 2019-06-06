# SOSTOOLS
A free MATLAB toolbox for formulating and solving sums of squares (SOS) optimization programs.

## Introduction
**SOSTOOLS** is a free MATLAB toolbox for formulating and solving sums of squares (SOS) optimization programs. SOSTOOLS can be used to specify and solve sum of squares polynomial problems using a very simple, flexible, and intuitive high-level notation. The SOS programs can be solved using [SeDuMi](http://sedumi.ie.lehigh.edu/), [SDPT3](http://www.math.nus.edu.sg/~mattohkc/sdpt3.html), [CSDP](https://projects.coin-or.org/Csdp/), [SDPNAL](http://www.math.nus.edu.sg/~mattohkc/SDPNAL.html), [SDPNAL+](http://www.math.nus.edu.sg/~mattohkc/SDPNALplus.html), [CDCS](https://github.com/oxfordcontrol/CDCS) and [SDPA](http://sdpa.sourceforge.net/). All these are well-known semidefinite programming solvers, with SOSTOOLS handling internally all the necessary reformulations and data conversion.

## What is a "sum of squares optimization program"? Why would I want such a thing?

A *sum of squares (SOS) program*, in the simplest case, has the form:
```
minimize: c_1 * u_1 + ... + c_n * u_n

subject to constraints: 

P_i(x) := A_i0(x) + A_i1(x) * u_1 + ... + A_in(x) * u_n

are sums of squares of polynomials (for i=1..n).
```

Here, the `A_ij(x)` are multivariate polynomials, and the decision variables `u_i` are scalars. This is a convex optimization problem, since the objective function is linear and the set of feasible `u_i` is convex.

While this looks quite nice, perhaps you are actually interested in more concrete problems such as:

* Constrained or unconstrained optimization of polynomial functions.
* Mixed continuous-discrete optimization.
* Finding Lyapunov or Bendixson-Dulac functions for nonlinear dynamical systems (with polynomial vector fields).
* Deciding copositivity of a matrix.
* Inequalities in probability theory.
* Distinguishing separable from entangled states in quantum systems.

Or, more generally, problems that deal with basic semialgebraic sets (sets defined by polynomial equalities and inequalities).

## Distribution and release information

SOSTOOLS is freely available under the GNU public license v3.0.

Archived releases of SOSTOOLS prior to v3.03 can be bounds at the following sites:

* CDS at Caltech: http://www.cds.caltech.edu/sostools 

* LIDS at MIT: http://www.mit.edu/~parrilo/sostools

* Control Group at Oxford: http://www.eng.ox.ac.uk/control/sostools

## System requirements

To install and run SOSTOOLS, you need:

* [MATLAB](http://www.mathworks.com/) version 6.0 or later.
* MATLAB Symbolic Math Toolbox version 2.1.2 (optional) for SOSTOOLS versions 2.05 and earlier, or the current version of the MATLAB Symbolic Math Toolbox for SOSTOOLS version 3.00 and later.
* An SDP solver, either [SeDuMi](http://sedumi.ie.lehigh.edu/), [SDPT3](http://www.math.nus.edu.sg/~mattohkc/sdpt3.html), [CSDP](https://projects.coin-or.org/Csdp/), [SDPNAL](http://www.math.nus.edu.sg/~mattohkc/SDPNAL.html), [SDPNAL+](http://www.math.nus.edu.sg/~mattohkc/SDPNALplus.html), [CDCS](https://github.com/oxfordcontrol/CDCS) and [SDPA](http://sdpa.sourceforge.net/). These solvers and their documentation can be downloaded for free. For information on how to install them, you are referred to their installation instructions.
* SOSTOOLS can easily be run on Windows or MAC OSX machines. It utilizes MATLAB sparse matrix representation for good performance and to reduce the amount of memory needed.

Detailed installation instructions are available in the SOSTOOLS user's guide [here](docs/sostools.pdf).

## Authors

The software has been written and is maintained by:

* [Antonis Papachristodoulou](http://sysos.eng.ox.ac.uk/control/sysos/index.php/User:Antonis)
* [James Anderson](http://sysos.eng.ox.ac.uk/control/sysos/index.php/User:James_Anderson)
* [Giorgio Valmorbida](http://www.l2s.centralesupelec.fr/perso/giorgio.valmorbida)
* [Stephen Prajna](http://www.cds.caltech.edu/~prajna/)
* [Peter Seiler](http://www.aem.umn.edu/people/faculty/bio/seiler.shtml)
* [Pablo A. Parrilo](http://www.mit.edu/~parrilo)

## References
For a detailed explanation of the theory and applications of sums of squares programming, as well as references to related work, please see:

* *Structured Semidefinite Programs and Semialgebraic Geometry Methods in Robustness and Optimization*.
California Institute of Technology, Pasadena, CA, May 2000. [Abstract](http://www.mit.edu/~parrilo/pubs/files/Thesis_abstract.html), [Pdf](http://www.mit.edu/~parrilo/pubs/files/thesis.pdf). 

* *Semidefinite programming relaxations for semialgebraic problems*.
P. A. Parrilo.  [Abstract](http://www.mit.edu/~parrilo/pubs/files/SDPrelax_abstract.html) [ps](http://www.mit.edu/~parrilo/pubs/files/SDPrelaxations.ps), [gz](http://www.mit.edu/~parrilo/pubs/files/SDPrelaxations.ps.gz).

* *Minimizing polynomial functions*
P. A. Parrilo, B. Sturmfels. [ArXiV](http://www.arxiv.org/abs/math.OC/0103170).


For more references please see http://hot.caltech.edu/math.html but also the authors' websites.
 

## Feedback
For comments, bug reports, encouragement, suggestions, complaints, etc., please send email to: sostools@cds.caltech.edu.

If you use SOSTOOLS for research purposes, we'd be happy to hear about it and mention it in the reference guide. Please drop us a line, to sostools@cds.caltech.edu.

## Citing 

Please use the following when citing SOSTOOLS:

```
@manual{sostools,
author = {A. Papachristodoulou, J. Anderson, G. Valmorbida, S. Prajna, P. Seiler and P. A. Parrilo},
title = {{SOSTOOLS}: Sum of squares optimization toolbox for {MATLAB}},
note = {Available from \texttt{http://www.eng.ox.ac.uk/control/sostools}, \texttt{http://www.cds.caltech.edu/sostools} and \texttt{http://www.mit.edu/\~{}parrilo/sostools}},
year = {2013},
address = {\texttt{http://arxiv.org/abs/1310.4716}},
}
```
## Related Patches and Add-Ons

The following are some patches and add-ons written by other people:
* [INTSOSTOOLS](https://github.com/gvalmorbida/INTSOSTOOLS) for formulating and solving optimization problems subject to one-dimensional integral inequalities.
* [frlib](https://github.com/frankpermenter/frlib) for a pre-processing facial reduction step.
