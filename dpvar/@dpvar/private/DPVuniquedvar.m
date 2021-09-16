function [newC,uniquedvarname] = DPVuniquedvar(C,dvarname)
% [newC,uniquedvarname] = DPVuniquedvar(C,dvarname) removes
% repeated decision variables from dvarname, and adds the associated
% rows in C to obtain a compressed C matrix.
%
% INPUTS:
% C: a coefficient matrix associated with some dpvar object, #rows of C
% should be exactly divisible by length(dvarname)+1
% dvarname: a cell of decision variable names, associated with the same dpvar
% 
% OUTPUTS:
% uniquedvarname: a cell with the same dvarnames as the input, but no repetitions
% newC: associated coefficient matrix with rows corresponding to the same
% dvarname added
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
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% 07/01/2021: DJ, MP, SS, Copied from polynomial PVuniquevar by PJS
%               Adjusted for use on dvarnames of dpvars


% Find repeated variables:
%  sortidx:   index such that varsort = dvarname(sortidx)
%  repeatidx: nonzero indices indicate repeats
%  rvec:      map entries of varsort to numbers 1,..,nud
[varsort,sortidx] = sortrows_integerTable(char(dvarname));
repeatidx = all(varsort(2:end,:) == varsort(1:end-1,:),2);
rvec = cumsum([1; ~repeatidx]);

% Create unique set of decision variables
uniqueidx = sortidx;
uniqueidx( repeatidx ) = [];
uniquedvarname = dvarname(uniqueidx);

% Form coefficient matrix in terms of new variable set
%  (Use sparse multiply, summat*C, to sum the rows)
% For dpvars, each row in C is associated with a dvar
% or a constant term, and a row in the final (matrix-valued) object.
% We have to adjust the indices accordingly.
nd = length(dvarname)+1;        % note +1 to account for constant term
nud = length(uniquedvarname)+1;
nrows = size(C,1);
nrows_new = nud * nrows/nd;
rvec_full = reshape([1;rvec+1] + (0:nud:nrows_new-nud),nrows,1);    % indices of rows in output C; note +1 to account for constant term
sortidx_full = reshape([1;sortidx+1] + (0:nd:nrows-nd),nrows,1);    % indices of rows in input C
summat = sparse(rvec_full,sortidx_full,1,nrows_new,nrows,nrows);
newC = summat * C;

end