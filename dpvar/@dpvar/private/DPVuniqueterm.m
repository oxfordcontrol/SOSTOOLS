function [newC,uniquedeg] = DPVuniqueterm(C,degmat)
% [newC,uniquedvarname] = DPVuniqueterm(C,degmat) removes
% repeated monomial terms from degmat, and combines the associated
% columns in C in accordance with dpvar structure
%
% INPUTS:
% C: a coefficient matrix associated with some dpvar object, # columns of C
% must be divisible by # rows of degmat
% degmat: a sparse array of monomial term degrees, associated with the same dpvar
% 
% OUTPUTS:
% uniquedeg: a new array with the same rows as degmat, but no repetitions
% newC: associated coefficient matrix with columns corresponding to the same
% monomial term added
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
% 07/01/2021: DJ, MP, SS Copied from polynomial PVuniqueterm by PJS
%               Adjusted for use on dpvars

% Find repeated monomials
% Sort done by total monomial degree first and then lexicographic
%   sortidx: index such that degsort = degmat(sortidx)
%   repeatidx: nonzero indices indicate repeats
%   rvec:  map rows of degsort to numbers 1,..,nut
M = [sum(degmat,2) degmat];
[degsort,sortidx] = sortrows_integerTable(M);
degsort = degsort(:,2:end);
repeatidx = ~any(degsort(2:end,:) ~= degsort(1:end-1,:),2);
rvec = cumsum([1; ~repeatidx]);

% Create unique set of monomials
uniqueidx = sortidx;
uniqueidx( repeatidx ) = [];
uniquedeg = degmat(uniqueidx,:);

% Sum repeated terms
%  (Use sparse multiply, C*summat, to sum the columns)
nt = size(degmat,1);
nut = size(uniquedeg,1);
% DJ 07/01/21: For dpvars, each column in C is associated with a row in
% degmat and column in the final (matrix-valued) object. We have to adjust
% rvec and sortidx accordingly.
ncols = size(C,2);
ncols_new = nut * ncols/nt;
rvec_full = reshape(rvec + (0:nut:ncols_new-nut),ncols,1);  % columns in output C
sortidx_full = reshape(sortidx + (0:nt:ncols-nt),ncols,1);  % columns in input C
summat = sparse(sortidx_full,rvec_full,1,ncols,ncols_new,ncols);
newC = C*summat;

end