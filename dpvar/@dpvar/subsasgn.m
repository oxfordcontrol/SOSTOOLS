function a = subsasgn(a,L,RHS)
% a = subsasgn(a,L,RHS) assigns a value RHS to property L of dpvar a
%
% INPUT
% a: dpvar object
% L: struct of variable size
%       - Each element should have two fields:
%           > type: either '.', or '()'
%           > subs: either a char, a cell, an array, or the value ':'
%       - If L(1).type='.', L(1).subs should be one of the fieldnames 
%         'matdim', 'degmat', 'C', 'varname', or 'dvarname' of
%         the dpvar object.
%       - If L(i).type='()', L(i).subs should be an array of indices
%         indicating desired elements of the object to be adjusted. If
%         ref(i).subs=':', all elements will be adjusted To specify row
%         and column indices, use L(i).subs{1} for row indices, and
%         L(i).subs{2} for column indices.
% RHS: dpvar, polynomial, array, cell or char
%       - Object class of RHS should match object class of subsref(a,L)
%         (i.e., if subsref(a,L) = a.dvarname, RHS should be a cell of
%         chars specifying decision variable names).
%         Exception to this is that elements, rows or columns of the dpvar
%         may also be specified by a polynomial or double, instead of a 
%         dpvar.
%       - Size of RHS should always match size of subsref(a,L).
% 
% OUTPUT
% a: dpvar object
%       - Object is adjusted so that subsref(a,L) = RHS
%       - If L(1).type='.', a will be such that field L(1).subs of the 
%         dpvar has a value RHS
%       - If L(1).type='()', a will be such that Dpsub will be a dpvar object containing only
%         the desired rows and columns of Dp as specified by ref(1).subs.
%
% EXAMPLE:
% Let Dp be a dpvar specifying a polynomial
% polynomial(Dp) = [d1, x; d2*y, x*y];
% Then:
%   - Typing Dp.dvarname = {'d3';'d4'} will call Dp = subsasgn(Dp,L,RHS)                     
%     with L.type = '.', L.subs = 'dvarname', and RHS = {'d3';'d4'},
%     so that polynomial(Dp) = [d3, x; d4*y, x*y];
%
%   - Typing Dp.dvarname{2} = 'd4' will call Dp = subsasgn(Dp,L,RHS)
%     with L(1).type = '.', L(1).subs = 'dvarname', 
%          L(2).type = '{}', L(2).subs = 2, and RHS = 'd4',
%     so that polynomial(Dp) = [d1, x; d4*y, x*y];
%
%   - Typing Dp(1,2) = 3 will call Dp = subsasgn(Dp,L,RHS)
%     with L.type = '()', L.subs{1} = 1, L.subs{2} = 2, and RHS = 3,
%     so that polynomial(Dp) = [d1, x; 3, x*y];
%     Note that RHS = 3 can be either a double, polynomial, or dpvar class
%     object.
%
%   - Typing Dp(:,2) = [3;d2] will call Dp = subsasgn(Dp,L,RHS)
%     with L.type = '()', L.subs{1} = ':', L.subs{2} = 2, and RHS = [3;x],
%     so that polynomial(Dp) = [d1, x; 3, d2];
%     Note that RHS = [3;d2] MUST be a dpvar object in this case
% 
% NOTES:
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding, DJ, MP, SS - 09/27/2021
% added correction to allow dynamic extension of dpvar size, SS - 8/10/2021
% 02/21/2022 - DJ: Update for case '.' to avoid use of "set".
% 10/31/2024 - DJ: Add support for linear indexing;
% 01/28/2026 - DJ: Replace sparse indexing with manually rebuilding matrix;

switch L(1).type
    case '.'
    % Assign a value RHS to one of the fields (e.g. a.matdim) of dpvar a
        if ~ismember(L(1).subs,{'dvarname','varname','C','degmat','matdim'})
            error(['Input "',L(1).subs,'" is not an assignable field of dpvar objects'])
        end
        if length(L) == 1
            temp = RHS;
        else
            % Peform all subsasgn but L(1)
            temp = subsref(a,L(1));
            temp = subsasgn(temp,L(2:end),RHS);
        end
        a.(L(1).subs) = temp;
    case '()'
    % Assign a value to certain rows/columns of dpvar a
        if length(L)~=1
            % disp("subassgn of the type 'dpvar(a,b).prop = RHS;' is currently unsupported");
            % Peform all subsasgn but L(1)
            temp = subsref(a,L(1));
            temp = subsasgn(temp,L(2:end),RHS);
            % Perform the subsasgn L(1)
            a = subsasgn(a,L(1),temp);
            return
        end
        
        if ~ isa(RHS,'dpvar') && ~isa(RHS,'polynomial') && ~isa(RHS,'double')
            error('New values of elements must be specified as object of type ''dpvar'', ''polynomial'', or ''double''.')
        end

        % % First, process the indices specifying which elements of a to
        % % change.
        if isscalar(L(1).subs) 
            % Assume linear indices to be specified.
            use_lindx = true;
            indl = L(1).subs{1};

            % Check that the indices make sense.
            if strcmp(indl,':')
                indl = 1:numel(a);
            elseif ~isnumeric(indl) || any(indl<=0)
                error('Indices must be strictly positive integers.')
            end
            
            % Allow indices outside of domain for vector valued setting.
            if any(indl>numel(a)) && any(a.matdim<=1)
                if a.matdim(2)==1 && a.matdim(1)~=1
                    % Grow column vector
                    a.matdim(1) = max(indl);
                else
                    % Grow row vector
                    a.matdim = [1,max(indl)];
                end
                % expand coefficient matrix with zeros to match new size. 
                tmp = a.C;
                [i,j,v] = find(tmp);
                a.C = sparse(i,j,v,a.matdim(1)*(length(a.dvarname)+1),a.matdim(2)*size(a.degmat,1),nnz(a.C));
            elseif any(indl>numel(a))
                error('Indices cannot exceed dimensions of the object.')
            end
        elseif length(L(1).subs)==2
            % Assume row and column indices to be specified.
            use_lindx = false;
            indr = L(1).subs{1}(:)';
            indc = L(1).subs{2}(:)';

            % Convert ':' to actual indices
            if strcmp(indr,':')
                indr = 1:a.matdim(1);
            elseif ~isnumeric(indr) || any(indr<=0)
                error('Indices must be strictly positive integers.')
            end
            if strcmp(indc,':')
                indc = 1:a.matdim(2);
            elseif ~isnumeric(indc) || any(indc<=0)
                error('Indices must be strictly positive integers.')
            end
            
            % Allow index exceeding matrix dimension by augmenting with 0.
            if ~isempty(RHS)
                if any(a.matdim==0)
                    a.matdim = [max(indr),max(indc)];
                end
                if any(indr>a.matdim(1))
                    a.matdim(1) = max(indr);
                end
                if any(indc>a.matdim(2))
                    a.matdim(2) = max(indc);
                end
                % expand coefficient matrix and reassign values 
                tmp = a.C;
                [i,j,v] = find(tmp);
                a.C = sparse(i,j,v,a.matdim(1)*(length(a.dvarname)+1),a.matdim(2)*size(a.degmat,1),nnz(a.C));
    
                % Convert row and column indices to linear ones.
                indr_f = repmat(indr,[1,length(indc)]);
                indc_f = kron(indc,ones(1,length(indr)));
                indl = sub2ind(size(a),indr_f,indc_f);
            end
        else
            error('Assigning values of elements of dpvar objects is only supported using linear indices or row and column indices.')
        end

        % Check that the dimensions of the indices and RHS match.
        if ~isempty(RHS)
            if numel(RHS)==1
                RHS = RHS*ones(1,numel(indl));
            elseif numel(RHS)~=numel(indl)
                error('Number of elements of which to change value should match number of new values.')
            end
            indl = indl(:)';
        end
         
        
        % % Next, process the RHS.
        % If RHS is empty, we assume the user wants to remove the specified
        % rows or columns from the object.
        if isempty(RHS)
            if use_lindx
                % Removing elements of "a" specified by linear indices
                % is supported only if "a" is a vector.
                if a.matdim(1)==1
                    indr = 1;       indc = indl;
                elseif a.matdim(2)==1
                    indr = indl;    indc = 1;
                else
                    error('Excluding elements of ''dpvar'' object using linear indexing is not supported.')
                end
            end
            a = col_row_remove(a,indr,indc);
            return   
        end
        
        % Otherwise, make sure "RHS" is of type 'dpvar', and that the
        % dvarnames, varnames, and degmats of "a" and "RHS" match.
        if isa(RHS,'dpvar') % merge bases
            [a,RHS] = common_basis(a,RHS);
        elseif isa(RHS,'polynomial') % convert to dpvar, then merge bases
            RHS = poly2dpvar(RHS);
            [a,RHS] = common_basis(a,RHS);
        elseif isa(RHS,'double') % convert to dpvar, no change needed for dvarname, varname or degmat
            RHS = dpvar(RHS);
            [a,RHS] = common_basis(a,RHS);
        else
            error('New values of elements must be specified as object of type ''dpvar''.')
        end
        
        % % Finally, perform the actual subsasgn.
        % To adjust the elements of "a", we need to adjust the coefficient
        % matrix a.C. For this, identify corresponding locations of 
        % coefficients in "a" that need to be changed.
        if use_lindx
            idxl_C = getCindices(a,indl);
        else
            idxl_C = getCindices(a,indr_f,indc_f);
        end
        % Also establish indices of elements in RHS.C associated
        % to each of the elements of RHS.
        idxl_C_RHS = getCindices(RHS,(1:numel(RHS)));
        % All that remains is to replace the coefficients of "a" with those
        % of "RHS".
        % We do this manually, extracting the nonzero elements from a.C,
        % removing elements at "idxl_C", and replacing them with the
        % elements from RHS
        %a.C(idxl_C(:)) = RHS.C(idxl_C_RHS);     % slow implementation...
        [r_idcs,c_idcs,vals] = find(a.C);
        l_idcs = sub2ind(size(a.C),r_idcs,c_idcs);
        l_retain = ~ismember(l_idcs,idxl_C(:));
        l_full = [l_idcs(l_retain); idxl_C(:)];
        val_full = [vals(l_retain); RHS.C(idxl_C_RHS(:))];
        [r_full,c_full] = ind2sub(size(a.C),l_full);
        a.C = sparse(r_full,c_full,val_full,size(a.C,1),size(a.C,2));       % DJ, 01/28/2026

    case '{}'
        error('{}- like subsassign is not supported for dpvar objects.');
end

% Recheck if all properties have correct type/dimension
[logval, errmsg] = dpvarconstructmsg(a);
if logval~=0
    a.chkval=0;
    error(errmsg);
else
    a.chkval=1;
end

end



function [indl_C,indc] = getCindices(Dp,ind1,ind2)
% For a given dpvar object Dp, and a 1xn vector of linear indices ind1,
% this returns a mxn array indl_C specifying for each of the n elements
% Dp(indl(i)) which linear indices indl_C(:,i) in the coefficient matrix
% Dp.C correspond to this particular element.
% Alternatively, indices can be specified as row indices ind1, and column
% indices ind2.

% Check number of decision variables and monomials.
Zdlen = length(Dp.dvarname) + 1;
Zplen = size(Dp.degmat,1);
% Determine row and column index in Dp associated to each linear index.
if nargin==2
    [indr,indc] = ind2sub(size(Dp),ind1);
else
    indr = ind1;    indc = ind2;
end
% For each element, determine which rows and columns in Dp.C correspond to
% this element.
indr_C = Zdlen*(indr-1) +(1:Zdlen)';
indc_C = Zplen*(indc-1) +(1:Zplen)';
% Extend arrays to account for fact that each row may appear in multiple
% columns, and each column may appear in multiple rows.
indr_C = repmat(indr_C,[Zplen,1]);
indc_C = kron(indc_C,ones(Zdlen,1));
% Convert back to linear indices;
indl_C = sub2ind(size(Dp.C),indr_C,indc_C);

end



function a = col_row_remove(a,indr,indc)
% Remove columns and rows of "a" specified by indr and indc.
% Note that either "indr" must specify full columns of "a", or "indc" must
% specify full rows of "a", as removing just a single element from a
% matrix-valued object is not supported.

% Check which rows and columns the user wants to retain.
rfull = 1:a.matdim(1);      rretain = setdiff(rfull,indr);
cfull = 1:a.matdim(2);      cretain = setdiff(cfull,indc);

if isempty(rretain) && isempty(cretain)
    % Set all elements of a equal to empty.
    a = dpvar([]);
elseif isempty(rretain)
    % Set columns specified by indc to empty.    
    Ctemp = a.C;
    nz = size(a.degmat,1);
    cCretain = cell2mat(cellfun(@(x) (x-1)*nz+1:x*nz,num2cell(cretain),'uni',0));
    a.C = Ctemp(:,cCretain);
    a.matdim = [a.matdim(1),length(cretain)];
    a = combine(a);
elseif isempty(cretain)
    % Set rows specified by indr to empty.
    Ctemp = a.C;
    nd = length(a.dvarname)+1;
    rCretain = cell2mat(cellfun(@(x) (x-1)*nd+1:x*nd,num2cell(rretain),'uni',0));
    a.C = Ctemp(rCretain,:);
    a.matdim = [length(rretain),a.matdim(2)];
    a = combine(a);
else
    error('Null assignment must be specified for entire row or column.')
end

end