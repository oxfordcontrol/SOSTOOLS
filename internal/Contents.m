% SOSTOOLS -- Sum of Squares Toolbox, internal files
% Version 4.00 
% 9 September 2021.
%
% General functions used internally by SOSTOOLS
% CONVERTTOCHAR     --- Convert symbolic object to char.
% FINDCOMMONZ       --- Combine monomial vectors.
% LRSCALE           --- Applies left and right scaling matrices.
% SORTNOREPEAT      --- Exclude rows from a matrix appearing in another matrix.
% SPARSEMULTIPART   --- Find elements in set of monomials that lie inside
%                       convex hull in set of bipartite monomials.
%
% - Functions introduced for dpvar implementation
% FASTINT2STR       --- Construct coefficient names from integer array.
% SOSPOLYVAR_MAT    --- Declare a new polynomial matrix variable in an SOSP.
%
% - Functions used for internal processing
% FRLIB_POST        --- Displays the reductions obtained with frlib.
% FRLIB_PRE         --- Calls frlib for reduction of an SOSP.
% GCDS              --- Computes greatest common denominator of a list of
%                       numbers.
% INCONVHULL        --- Find elements in set of monomials that lie inside
%                       convex hull of another set.
% LCMS              --- Computes least common multiple of a list of numbers.
% MOSEKSOL2SEDUMISOL--- Convert Mosek output to Sedumi output.
% MYSVEC            --- Convert symmetric matrix to vector.
% PROJ3             --- Perform orthogonal projection of point onto hyperplane.
% SEDUMI2MOSEK      --- Convert Sedumi input to Mosek input.
% USECONVHULL       --- Find defining inequalities of convex hull of set of 
%                       monomials.
% VREP2HREP         --- Compute H-representation of convex hull.
%
% - Functions used for building the SOS program (SOSP) strucuture
% GETCONSTRAINT     --- Find constraint for sum of squares decomposition.
% GETEQUATION       --- Convert symbolic constraint to SOSP format.
% SOSCONSTR         --- Add a new constraint to an SOSP.
% SOSVAR            --- Declare a new variable in an SOSP.
%
% - Functions for faster operations exploiting structure and sparsity
% MAKESPARSE                --- Convert matrix to sparse format
% SORTROWS_INTEGERTABLE     --- Sort rows of array of integers
% SPANTIBLKDIAG             --- Anti block diagonal concatenation of sparse
%                               matrices
% SPBLKDIAG                 --- Block diagonal concatenation of sparse
%                               matrices
% SPREPMAT                  --- Repmat for sparse matrices
%
% - Functions to handle symbolic expressions
% GETDEGREES        --- Get degrees of a monomial
% GETPOLYSYM        --- Convert symbolic polynomial to SOSP format
% MYSYMPOWER        --- Take power of symbolic object
% MYSYMSUBS         --- Symbolic substitution
% SYM2CHARTABLE     --- Convert symbolic vector to string table

% This file is part of SOSTOOLS - Sum of Squares Toolbox ver 4.00.
%
% Copyright (C)2002, 2004, 2013, 2016, 2018, 2021  A. Papachristodoulou (1), J. Anderson (1),
%                                      G. Valmorbida (2), S. Prajna (3), 
%                                      P. Seiler (4), P. A. Parrilo (5)
%                                      M. Peet (6), D. Jagt (6)
% (1) Department of Engineering Science, University of Oxford, Oxford, U.K.
% (2) Laboratoire de Signaux et Systmes, CentraleSupelec, Gif sur Yvette,
%     91192, France
% (3) Control and Dynamical Systems - California Institute of Technology,
%     Pasadena, CA 91125, USA.
% (4) Aerospace and Engineering Mechanics Department, University of
%     Minnesota, Minneapolis, MN 55455-0153, USA.
% (5) Laboratory for Information and Decision Systems, M.I.T.,
%     Massachusetts, MA 02139-4307
% (6) Cybernetic Systems and Controls Laboratory, Arizona State University,
%     Tempe, AZ 85287-6106, USA.
%
% Send bug reports and feedback to: sostools@cds.caltech.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.