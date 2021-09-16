function C = mtimes(E,F)
% C = mtimes(E,F) multiplies a dvpar object and a pvar polynomial object.
%
% INPUTS:
% E: A dpvar or pvar variable (only one object can be a dpvar), dimension
% [m n]
% F: A dpvar or pvar variable (only one object can be a dpvar), dimension [n
% p]
%
% OUTPUTS:
% C: a dpvar object, dimension [m p]
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
% Initial coding DJ, MP,SS - 07/18/2021

% If either object is double, convert to polynomial
if ~isa(E,'dpvar') && ~isa(E,'polynomial')
    E = polynomial(E);
end
if ~isa(F,'dpvar') && ~isa(F,'polynomial')
    F = polynomial(F);
end

% Start multiplication process, find which object is dpvar
if isa(E,'dpvar')&&isa(F,'dpvar') % 2 dpvars should not be multiplied
    error('Cannot multiply two dpvars');
    
elseif isa(F,'dpvar') % right object is dpvar --> multiply transposes instead
    C_transpose = F.'*E.';
    if all(size(C_transpose.C)==0) % correction for zero dimension dpvars
        C_transpose.C = sparse(size(C_transpose,1)*(length(C_transpose.dvarname)+1), size(C_transpose,2)*size(C_transpose.degmat,1));
    end
    C = C_transpose.';
    
else  % left object is dpvar
    Edegmat=E.degmat;       Fdegmat=F.degmat; % we assume there are no repeats in degmat...
    Epvarname = E.varname;  Fpvarname = F.varname;
    ntE=size(Edegmat,1);    ntF=size(Fdegmat,1);
    
    % synchronize the pvarnames
    [pvarnames_all,IE_p,~] = union(Epvarname,Fpvarname,'stable');   % returns indices in E
    [~,~,IF_p] = intersect(Fpvarname,pvarnames_all,'stable');       % returns indices in F
    np_all=length(pvarnames_all);
    
    % adjust the degmats accordingly
    Edegmat_new=spalloc(ntE,np_all,nnz(Edegmat));
    Fdegmat_new=spalloc(ntF,np_all,nnz(Fdegmat));
    Edegmat_new(:,IE_p)=Edegmat;
    Fdegmat_new(:,IF_p)=Fdegmat;
    
    % adjust E and F to use the new varnames and degmats
    F = polynomial(F.coef,Fdegmat_new,pvarnames_all,F.matdim);
    E.degmat=Edegmat_new;
    E.varname = vec(pvarnames_all); % NOTE: if E and F have only 1 varname, need vec to get appropriate dim
    
    % Distinguish three cases:
    if all(size(E)==[1,1]) % scalar dpvar multiplying a pvar
        C = E;
        C.matdim = size(F);
        
        % multiply coefficient matrices
        fcoefs = spalloc(size(F,1),size(F.coefficient,1)*size(F,2),nnz(F.coefficient));
        for i = 1:size(F,2) % maybe implement fancy reshape from pvar mtimes??
            tmp = F.coefficient(:,[(i-1)*size(F,1)+1:i*size(F,1)]);
            fcoefs(:,(i-1)*size(F.coefficient,1)+1:i*size(F.coefficient,1)) = tmp';
        end
%     Alternative implementations
%         fcoefs = [];
%         for i = 1:size(F,2) % maybe implement fancy reshape from pvar mtimes??
%             tmp = F.coefficient(:,[(i-1)*size(F,1)+1:i*size(F,1)]);
%             fcoefs = [fcoefs, tmp'];
%         end
%     Alternative implementation 2
%         fcoefs = Bshape(F); % this is faster only for large dimensions
        C.C = kron(fcoefs,E.C); % Construction of Cnew
        
        % multiply degmats
        %   this assumes both E and F have same pvars
        tmpdegmat = rowwisesum(C.degmat,F.degmat);
        
        % combine the product degmat
        [Adeg, degmat] = combine_degmat(tmpdegmat);
        
        % adjust coefficient matrix to combine rows with same monomials
        C.C = C.C*kron(eye(size(F,2)),Adeg);
        C.degmat = degmat;
        
    elseif all(size(F)==[1,1]) % matrix dpvar times scalar polynomial
        C = E;
        C.matdim = size(E);
        
        % multiply coefficient matrices
        fcoef = F.coefficient;
        ecoef = E.C;
        C.C = kron(ecoef,fcoef'); % Construction of Cnew
        
        % multiply degmats
        %   this assumes both E and F have same pvars
        tmpdegmat = rowwisesum(F.degmat,C.degmat);
        
        % combine the product degmat
        [Adeg, degmat] = combine_degmat(tmpdegmat);
        
        % adjust coefficient matrix to combine rows with same monomials
        C.C = C.C*kron(eye(size(E,2)),Adeg);
        C.degmat = degmat;
        
    else % dpvar matrix (E) times pvar matrix (F)
        if size(E,2)~=size(F,1)
            error('Inner dimensions dont match');
        end
        C = E;
        C.matdim = [E.matdim(1),size(F,2)];
        
        % multiply coefficient matrices
        fcoefs = Bshape(F); % this is faster only for large dimensions
        Cnew = C.C*kron(fcoefs,speye(size(E.degmat,1)));

        % multiply degmats
        %   this assumes both E and F have same pvars
        tmpdegmat = rowwisesum(C.degmat,F.degmat);
        
        % combine the product degmat
        [Adeg, degmat] = combine_degmat(tmpdegmat);
        
        % adjust coefficient matrix to combine rows with same monomials
        C.C = Cnew*kron(eye(size(F,2)),Adeg);
        C.degmat = degmat;
    end
end
end



function [A,d] = combine_degmat(dmat)
% Combine non-unqiue rows of dmat to build d, with A such that A*d=dmat
[d, ~, idxd] = unique(dmat,'rows');
A = sparse(1:length(idxd),idxd,ones(length(idxd),1),size(dmat,1), size(d,1));
% alternative way
% A = zeros(size(dmat,1), size(d,1));
% 
% for j=1:length(idxd)
%     A(j,idxd(j)) = 1;
% end
end
function [ C ] = rowwisesum( A, B ) 
% return the elementwise sum of two matrices, adds matrix A to every
% row of B 
C = kron(B,ones(size(A,1),1)) + repmat(A,size(B,1),1);
end