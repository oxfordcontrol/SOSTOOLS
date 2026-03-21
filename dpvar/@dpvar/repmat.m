function B = repmat(A,r1,r2)
% B = REPMAT(A,R1,R2) returns a the concatenation of R1 x R2 copies of A
% along the row and column directions
%
% INPUTS
% - A:  m x n dpvar object;
% - r1: scalar integer specifying the number of times A is to be repeated 
%       along the row direction;
% - r2: scalar integer specifying the number of times A is to be repeated
%       along the column direction;
%
% OUTPUTS
% - B:  m*r1 x n*r2 dpvar object representing kron(ones(nr,nc),A);
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2026  M. Peet, S. Shivakumar, D. Jagt
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
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% DJ, 03/20/2026: Initial coding;

% Check the inputs
if nargin<2
    error("Not enough input arguments.")
elseif nargin==2
    % Specify repetition factors as row vector
    repdim = r1;
elseif nargin==3
    repdim = [r1,r2];
else
    error("Replication is only supported along row and column directions.")
end
if any(size(repdim)~=[1,2])
    error("Replication factors must be specified as 1x2 vector of integers, or as 2 integer scalars.")
elseif any(~isnumeric(repdim)) || any(repdim~=round(abs(repdim)))
    error("Replication factors must be specified as nonnegative integers.")
end

% Split the dpvar object into individual components
matdim = A.matdim;
dvarname = A.dvarname;
varname = A.varname;
degmat = A.degmat;
C = A.C;

% Repeat the coefficients as desired
C_new = repmat(C,repdim);
matdim_new = matdim.*repdim;

% Return the repeated array
B = dpvar(C_new,degmat,varname,dvarname,matdim_new);  

end