function Z = monomials(vartable,d)
% MONOMIALS --- Construct a vector of monomials with 
%       prespecified degrees.
%
% Z = monomials(VARTABLE,DEGREE)
%
% Given a vector of independent variables VARTABLE (can be either
% symbolic or polynomial objects) and a vector of non-negative 
% integers DEGREE, this function constructs a column vector 
% containing all possible monomials of degree described in DEGREE.
%
% For example, monomials([x1,x2],[1:3]) with x1 and x2 
% being symbolic variables will return monomials in x1 and 
% x2 of degree 1, 2, and 3.
%

% Alternative calling convention:
%
%  Z = monomials(n,d)
%
% Given two integers n and d, this function constructs 
% a vector containing all possible monomials of degree d 
% in n variables.
%
% n is the number of variables
% d is the degree of monomials
% Z is the monomial vector
%

% This file is part of SOSTOOLS - Sum of Squares Toolbox ver 4.00.
%
% Copyright (C)2002, 2004, 2013, 2016, 2018, 2021  
%                                      A. Papachristodoulou (1), J. Anderson (1),
%                                      G. Valmorbida (2), S. Prajna (3), 
%                                      P. Seiler (4), P. A. Parrilo (5),
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



% 12/08/01 - SP
% 03/01/02 - SP -- Symbolic monomial


i = find(d<0);
if ~isempty(i)
    error('Negative degree.');
end;

Z = [];
for i = d
    Z = [Z; oldconstructZ(vartable,i)];
end;

% ==================================================================
function Z = oldconstructZ(vartable,d)
% Old constructZ

if isnumeric(vartable)
    n = vartable;
else
    n = length(vartable);
end;

ZZ = sparse(1,n);
for i = 1:n
    ss = size(ZZ,1);
    ZZ = sprepmat(ZZ,d+1,1);
    for j = 0:d
        ZZ(ss*j+1:ss*j+ss,i) = j;
    end;
    idx = find(sum(ZZ,2) <= d);   % Throw away invalid monomials
    ZZ = ZZ(idx,:);
end;
idx = find(sum(ZZ,2) == d);
Z = ZZ(idx,:);

if isa(vartable,'sym')
    vartable  = reshape(vartable,1,length(vartable));%makes vartable a row vector
    Z = prod(kron(vartable,ones(size(Z,1),1)).^Z,2);
elseif isa(vartable,'polynomial')
  coefftemp = speye(size(Z,1));
  Z = polynomial(coefftemp,Z,vartable.varname,[size(Z,1),1]);
end;


        
