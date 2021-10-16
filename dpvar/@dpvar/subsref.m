function Dpsub=subsref(Dp,ref)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function allows subindexing of the components of dpvar or extracts
% block matrices of dpvar. This function also allows accessing the
% properties of the dpvar class object using dot indexing
% 
% INPUT
% Dp: dpvar object
% ref: struct of variable size
%       - Each element should have two fields:
%           > type: either '.', or '()'
%           > subs: either a char, a cell, an array, or the value ':'
%       - If ref(1).type='.', ref(1).subs should be one of the fieldnames 
%         'matdim', 'degmat', 'C', 'varname', 'dvarname', or 'chkval' of
%         the dpvar object.
%       - If ref(i).type='()', ref(i).subs should be an array of indices
%         indicating desired elements to extract from the object. If
%         ref(i).subs=':', all elements will be extracted. To specify row
%         and column indices, use ref(i).subs{1} for row indices, and
%         ref(i).subs{2} for column indices.
% 
% OUTPUT
% Dpsub: cell, array, or dpvar object
%       - If ref(1).type='.', Dpsub will be the desired field of the dpvar
%         object (or elements of this field as specified by ref(2:end))
%       - If ref(1).type='()', Dpsub will be a dpvar object containing only
%         the desired rows and columns of Dp as specified by ref(1).subs.
%
% EXAMPLE:
% Let Dp be a dpvar specifying a polynomial
% polynomial(Dp) = [d1, x; d2*y, x*y];
% Then:
%   - Dp.dvarname = subsref(Dp,ref) = {'d1';'d2'}
%           where ref.type = '.' and ref.subs = 'dvarname'
%   - Dp.dvarname{2} = subsref(Dp,ref) = 'd2'
%           where ref(1).type = '.' and ref(1).subs = 'dvarname' and
%                 ref(2).type = '{}' and ref(2).subs = 2
%   - Dp.matdim(1) = subsref(Dp,ref) = 2
%           where ref(1).type = '.' and ref(1).subs = 'matdim' and
%                 ref(2).type = '()' and ref(2).subs = 1
%   - Dp(1,:) = subsref(Dp,ref) = [d1, x]
%           where ref.type = '()' and ref.subs{1} = 1 and ref.subs{2} = ':'
%   - Dp(4) = subsref(Dp,ref) = [x*y]
%           where ref.type = '()' and ref.subs = 4 (note that linear
%           indexing is used)
%  
% NOTES:
% - In general, there is no need to call subsref explicitly. Just type,
%   e.g. Dpsub = Dp(4,5), and Matlab will call subsref with ref.type =
%   '()', ref.subs{1} = 4, and ref.subs{2} = 5 for you.
%
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
% Initial coding DJ, MP, SS - 07/14/2021
% bug fix correction for linear indexing Dp(:), index exceeds dimension -
% SS 8/10/2021
% bug fix correction for linear indexing Dp([1,2,3...]), incorrect
% dimension in sosineq, SS - 8/11/2021
% Bug fix for linear indexing, DJ - 10/15/2021

switch ref(1).type
    case '.'
        % Extract just one of the fields (e.g. Dp.matdim) of Dp
        Dpsub = get(Dp,ref);
    case '()'
        % Extract certain rows/columns of the matrix-valued Dp
        matdim = Dp.matdim;
        
        % Extract size of decision monomial and independent variable
        % monomial
        Zdlen = length(Dp.dvarname) + 1;
        Zplen = size(Dp.degmat,1);
        
        if length(ref(1).subs)==1
            % % Distinguish case of linear indexing
            
            if strcmp(ref(1).subs{1},':')
                ref(1).subs{1} = 1:matdim(1)*matdim(2);
            end
            
            % Error checking
            if any(ref(1).subs{1}<=0)
                error('Indices must be positive integers');
            elseif any(ref(1).subs{1}>matdim(1)*matdim(2))
                error('Indices exceed matrix dimensions');
            end
            
            % Output will be a column vector
            mdim = [length(ref(1).subs{1}),1];
            [indr,indc] = ind2sub(size(Dp),ref(1).subs{1});
            
            % Determine rows in coefficient matrix associated with the desired
            % rows of the actual object
            nr = length(indr);
            indr = kron(vec((indr-1)*Zdlen),ones(Zdlen,1)) + repmat((1:Zdlen)',nr,1);
            indr = repmat(indr,Zplen,1);
            
            % Determine columns in coefficient matrix associated with the desired
            % rows of the actual object
            indc = kron((vec(indc) - 1)*Zplen + (1:Zplen),ones(Zdlen,1));
            indc = indc(:);
            
            % Extract the appropriate elements from the coefficient matrix
            linidx = sub2ind(size(Dp.C),indr,indc);
            Csub = spalloc(mdim(1)*Zdlen,Zplen,nnz(Dp.C));
            Csub(:) = Dp.C(linidx);
            
            % Build the new dpvar with the adjusted coefficient matrix and dimension
            Dpsub = compress(dpvar(Csub, Dp.degmat, Dp.varname, Dp.dvarname,mdim));
            return
        
        elseif length(ref(1).subs)==2
            % % Consider case of two indices
            
            indr = ref(1).subs{1};
            indc = ref(1).subs{2};
            
            mdim = [0 0];
            % Determine rows in coefficient matrix associated with the desired
            % rows of the actual object
            if strcmp(indr,':')
                indr = 1:size(Dp.C,1);
                mdim(1) = matdim(1);
            elseif length(indr)>=2
                mdim(1) = length(indr);
                indr = kron(vec((indr-1)*Zdlen),ones(Zdlen,1)) + repmat((1:Zdlen)',mdim(1),1);
            elseif length(indr)==1
                indr = (indr(1)-1)*Zdlen+1:(indr(1))*Zdlen;
                mdim(1) = 1;
            else
                indr = [];
                mdim(1) = 0;
            end
            % Determine columns in coefficient matrix associated with the desired
            % rows of the actual object
            if strcmp(indc,':')
                indc = 1:size(Dp.C,2);
                mdim(2) = matdim(2);
            elseif length(indc)>=2
                mdim(2) = length(indc);
                indc = kron(vec((indc-1)*Zplen),ones(Zplen,1)) + repmat((1:Zplen)',mdim(2),1);
            elseif length(indc)==1
                indc = (indc(1)-1)*(Zplen)+1:(indc(1))*(Zplen);
                mdim(2) = 1;
            else
                indc = [];
                mdim(2) = 0;
            end
            
            % Error checking
            if any(indr<=0) || any(indc<=0)
                error('Indices must be positive integers');
            elseif any(indr>size(Dp.C,1))||any(indc>size(Dp.C,2))
                error('Indices exceed matrix dimensions');
            end
            
            % If the coefficient matrices are very sparse, extract only nonzero rows/columns
            if nnz(Dp.C) <= 10*length(Dp.dvarname)
                [Ci,Cj] = find(Dp.C);
                indxr_n = ismember(indr,Ci);
                indxc_n = ismember(indc,Cj);
                Csub = spalloc(length(indr),length(indc),sum(indxr_n)*sum(indxc_n));
                Csub(indxr_n,indxc_n) = Dp.C(indr(indxr_n),indc(indxc_n));
            else
                Csub = Dp.C;
                Csub = Csub(indr,indc);
            end
                        
            % Build the new dpvar with the adjusted coefficient matrix and dimension
            Dpsub = compress(dpvar(Csub, Dp.degmat, Dp.varname, Dp.dvarname,mdim));
            
        else
            % In case someone specifies more than two arguments            
            if all(cell2mat(ref(1).subs(3:end))==1)
                ref(1).subs = ref(1).subs(1:2);
                Dpsub = subsref(Dp,ref);
            else            
                error('At most two indices are supported for dpvar class objects; dpvars are 2D!')
            end
            
        end
end
