function [B] = Bshape(P)
% B = Bshape(B) converts a coefficient matrix in polynomial format to dpvar
% format
% [b11 b21 b31 b12 b22 b32 b13 b23 b33] ---> [b11' b12' b13'
%                                             b21' b22' b23'
%                                             b31' b32' b33']
% INPUTS:
% B: polynomial class object
% 
% OUTPUTS:
% B: reshaped coefficient matrix
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
% Initial coding DJ, MP, SS - 08/07/2021

m = P.matdim(1); n = P.matdim(2);
ntp=size(P.degmat,1);

B1=P.coeff';


if m*n==0 % correction for zero dimensional polynomial matrix
    B = zeros(m,n*ntp);
    return
end

% first convert to cell block matrix vector of size matdim(2)x1
Bc=mat2cell(B1,m*ones(n,1),ntp );
B=cell2mat(reshape(Bc,1,n));% reshape cells and return to matrix


% Alternative Implementation,
% B3=spalloc(m,n*ntp,nnz(B1));
% JB1=repmat([1:ntp],1,m*n); % scan each row of B1 then move to next column
% IB1=kron([1:m*n],ones(1,ntp)); % the next column
% 
% JB3=repmat(repmat([1:ntp],1,m),1,n) + kron([0:n-1],ones(1,ntp*m))*ntp; % place in same col for m times
% IB3=repmat(kron([1:m],ones(1,ntp)),1,n) ; 
% % repeat ith row nt times ,then move to next row. What you reach row m,
% % repeat whole process
% idxB1=sub2ind([m*n,ntp],IB1,JB1);
% idxB3=sub2ind([m,ntp*n],IB3,JB3);
% B3(idxB3)=B1(idxB1);

% Alternative implementation 2
% % transpose each element of the block matrix and convert back to a matrix
% if ~outtype
%     B = cell2mat(cellfun(@transpose, B, 'un', 0));
% else
%     B = cellfun(@transpose, B, 'un', 0);
% end
end