function [newdegmat,uniquevar] = DPVuniquepvar(degmat,varname)
% [newdegmat,uniquevar] = DPVuniquepvar(degmat,varname) removes
% superfluous independent variables from varname, as well as the associated
% columns in degmat, in accordance with dpvar structure
%
% INPUTS:
% degmat: an array describing monomials in the variables varname, associated with some dpvar object
% varname: a cell of independent variable names, associated with the same dpvar
% 
% OUTPUTS:
% uniquevar: a new cell with the same names as varname, but no repetitions
% newdegmat: associated degmat, with columns corresponding to the same
% varname added
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% S. Shivakumar at sshivak8@asu.edu, or D. Jagt at djagt@asu.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2021  M. Peet, S. Shivakumar, D. Jagt
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
% 07/01/2021: DJ Copied from polynomial PVuniquevar by PJS

% Find repeated variables:
%  sortidx: index such that varsort = varname(sortidx)
%  repeatidx: nonzero indices indicate repeats
%  rvec: map entries of varsort to numbers 1,..,nuv
[varsort,sortidx] = sortrows_integerTable(char(varname));
repeatidx = all(varsort(2:end,:) == varsort(1:end-1,:),2);
rvec = cumsum([1; ~repeatidx]);

% Create unique set of variables
uniqueidx = sortidx;
uniqueidx( repeatidx ) = [];
uniquevar = varname(uniqueidx);

% Form degree matrix in terms of new variable set
%  (Use sparse multiply, degmat*summat, to sum the columns)
nv = size(degmat,2);
nuv = length(uniquevar);
summat = spalloc(nv,nuv,nv);
idx = sub2ind([nv nuv],sortidx,rvec);
summat(idx) = 1;
newdegmat = degmat*summat;
    
end