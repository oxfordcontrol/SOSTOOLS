function Dt = transpose(D)
% Dt = transpose(D) performs matrix transpose of a dpvar object. If D.C is
% certain to have highly sparse structure, use transpose_sparse() for
% faster computation
% 
% INPUTS:
% D: dpvar object
% 
% OUTPUTS:
% Dt: dpvar object corresponding to the (real) matrix transpose of D, Dt = D.'
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% S. Shivakumar at sshivak8@asu.edu, or Declan Jagt at djagt@asu.edu

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
% Initial coding DJ, MP, SS  - 07/18/2021

% Check the density of the coefficient matrices, and optionally use
% sparse implementation
dnstyD = nnz(D.C)/(size(D.C,1)*size(D.C,2));
if dnstyD<=0.1
    Dt = DP_sparse_transpose(D);
    return
end

% Note that the degmat, varnames, and dvarnames of Dt are identical to
% those of D
Dt = D;
%     Dt.dvarname = D.dvarname;
%     Dt.varname = D.varname;
%     Dt.degmat = D.degmat;

% The matrix dimension, of course, need to be flipped: m <-> n
Dt.matdim = fliplr(D.matdim);


% % % Now we need to adjust the coefficient matrix % % %
% This matrix consists of m x n blocks, each associated with one
% element of the mxn matrix-valued polynomial:
% D.C = [C11, C12, ..., C1n;]
%       [C21, C22, ..., C2n;]
%       [ : ,  : ,  . ,  : ;]
%       [Cm1, Cm2, ..., Cmn ]
% In order to tranpose the matrix-valued polynomial, we need to transpose
% the block matrix, without transposing the actual blocks:
% Dt.C = [C11, C21, ..., Cm1;]
%        [C12, C22, ..., Cm2;]
%        [ : ,  : ,  . ,  : ;]
%        [C1n, C2n, ..., Cmn ]

% We start by transposing the full coefficient matrix
DtC = (D.C)';
% DtC = [C11', C21', ..., Cm1';]
%       [C12', C22', ..., Cm2';]
%       [ :  ,  :  ,  . ,  :  ;]
%       [C1n', C2n', ..., Cmn' ]

% Next, we convert the array DtC to a corresponding cell to separate the
% different blocks
DtC = mat2cell(DtC, size(D.degmat,1)*ones(Dt.matdim(1),1), (length(D.dvarname)+1)*ones(Dt.matdim(2),1));
% DtC = [{C11'}, {C21'}, ..., {Cm1'};]
%       [{C12'}, {C22'}, ..., {Cm2'};]
%       [  :  ,    :  ,   . ,   :   ;]
%       [{C1n'}, {C2n'}, ..., {Cmn'} ]

% Finally, we transpose each individual block element of Dtc, and convert 
% back to a matrix
DtC = cell2mat(cellfun(@transpose, DtC, 'un',0));
if all(size(DtC)==0) % correction to include dpvars with zero rows/cols
    DtC = sparse(size(Dt,1)*(length(Dt.dvarname)+1),size(Dt,2)*size(Dt.degmat,1));
end
% DtC = [C11, C21, ..., Cm1;]
%        [C12, C22, ..., Cm2;]
%        [ : ,  : ,  . ,  : ;]
%        [C1n, C2n, ..., Cmn ]

% With that, we can finish the tranposed dpvar
Dt.C = DtC;
end



function Dt = DP_sparse_transpose(D)
% Dt = DP_sparse_transpose(D) performs matrix transpose of a dpvar object 
% 
% INPUTS:
% D: dpvar object
%       - transpose_sparse performs better when the coefficient matrix D.C
%         is sparser. Density of less than 10% is used as (heuristic) bound
%         in "transpose.m"
% 
% OUTPUTS:
% Dt: dpvar object corresponding to the (real) matrix transpose of D, Dt = D.'
% 
% NOTES:
% - This implementation is an alternative to the original transpose
%   function, found to work better when D.C is (very) sparse
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding DJ, MP, SS - 06/11/2021

% Note that the degmat, varnames, and dvarnames of Dt are identical to
% those of D
Dt = D;

% The matrix dimension, of course, need to be flipped: m <-> n
Dt.matdim = fliplr(D.matdim);

% Extract the number of decision vars (+1 to account for a constant
% term), and number of monomial terms described by the degmat
ld = length(D.dvarname)+1;
lp = size(D.degmat,1);

% % % Now we need to adjust the coefficient matrix % % %
% To do so, we note that each coefficient D.C(i*m+k, j*n+l) is associated 
% with a decision variable k, a monomial l, and a particular row and column 
% (i,j) in the mxn matrix.
% All we have to do is swap the associated row and column index i<->j

% First, find non-zero elements and row/column locations of coefficients
[idxI, idxJ, val] = find(D.C);

% Then, find associated column Qj and row Qi in the matrix
Qj = ceil(idxJ/lp)-1;   
Qi = ceil(idxI/ld)-1;

% Also, find associated monomial index Rl, and decision variable index Rk
Rl = rem(idxJ,lp);  Rl(~Rl) = lp; %<- this sets zero values to lp, since remainder zero means exactly divisible
Rk = rem(idxI,ld);  Rk(~Rk) = ld;

% Now, we just switch the row and column indices, but keep associated dvars
% and monomials the same
idxI_t = Qj*ld+Rk; % associated dvar remains the same, column number becomes row number
idxJ_t = Qi*lp+Rl; % associated monomial remains the same, row number becomes column number

% Finally, build the new coefficient matrix with swapped rows/columns
DtC = sparse(idxI_t,idxJ_t,val,ld*Dt.matdim(1),lp*Dt.matdim(2));

% and finish the transposed dpvar
Dt.C = DtC;
end