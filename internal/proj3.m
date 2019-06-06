function [x,NN] = proj3(x0,A,b,N); 
%
% Orthogonal projection of the point x0 onto the subspace Ax=b 
%
% Here, all rational data
% A is assumed fat, and full row rank
% A*A' is diagonal
%
% All integer computations
% The solution is x/NN, the input is x0/N
%
% Should revise it at some point.

% This file is part of SOSTOOLS - Sum of Squares Toolbox ver 3.03.
%
% Copyright (C)2002, 2004, 2013, 2016, 2018  A. Papachristodoulou (1), J. Anderson (1),
%                                      G. Valmorbida (2), S. Prajna (3), 
%                                      P. Seiler (4), P. A. Parrilo (5)
% (1) Department of Engineering Science, University of Oxford, Oxford, U.K.
% (2) Laboratoire de Signaux et Systmes, CentraleSupelec, Gif sur Yvette,
%     91192, France
% (3) Control and Dynamical Systems - California Institute of Technology,
%     Pasadena, CA 91125, USA.
% (4) Aerospace and Engineering Mechanics Department, University of
%     Minnesota, Minneapolis, MN 55455-0153, USA.
% (5) Laboratory for Information and Decision Systems, M.I.T.,
%     Massachusetts, MA 02139-4307
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

%aa = (A*A');
aa = sum(A.*A,2);
N2 = lcms(unique(aa));
x = N2*x0 - A'*( (N2 ./aa) .* (A*x0-b*N)) ;

NN = N2*N;

% Make sure it's minimal
[dum1,dum2,nz]=find(unique(abs(x)));
cfact = gcds([nz ; NN]);
x = x/cfact ;
NN = NN/cfact ;
