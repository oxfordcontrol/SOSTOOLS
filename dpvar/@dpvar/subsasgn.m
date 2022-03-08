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
% - At this time, linear indexing is not supported, so that e.g. Dp(3)=b 
%   will set row 3 of Dp equal to b, rather than setting the third element 
%   equal to b.
% - At this time, removing elements of Dp by specifying e.g. Dp(3,:)=[] (to
%   remove the third row of the dpvar) is not supported
% - At this time, setting all elements of a row or column of Dp equal to
%   the same value using scalar RHS (e.g. Dp(3,:) = 1) is not supported. 
%   Instead, make sure dimensions of RHS match those of subsref(a,L).
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
        %a = set(a,L(1).subs,temp); % <-- We don't have "set" for dpvar...
        a.(L(1).subs) = temp;
    case '()'
    % Assign a value to certain rows/columns of dpvar a
        if length(L)==1
            % Four '()'-subsasgn cases
            % In each case, make sure the dvarnames, varnames, and degmats
            % of "a" and "RHS" match using DPcommon_basis
            if isempty(a)
                a = dpvar();
            end
            if isempty(RHS)
                a = col_row_remove(a,L);
                return            
            elseif isa(RHS,'dpvar') % merge bases
                [a,RHS] = common_basis(a,RHS);
            elseif isa(RHS,'polynomial') % convert to dpvar, then merge bases
                RHS = poly2dpvar(RHS);
                [a,RHS] = common_basis(a,RHS);
            elseif isa(RHS,'double') % convert to dpvar, no change needed for dvarname, varname or degmat
                RHS = dpvar(RHS);
                [a,RHS] = common_basis(a,RHS);
            end
            
            % Next, determine the rows and columns of the matrix that need to be changed
            if length(L(1).subs)==1 
                warning("dpvar subassgn() does not support linear indexing. Assuming the first index corresponds to row locations and all columns");
                indr = L(1).subs{1};
                indc = ':';
            elseif length(L(1).subs)==2
                indr = L(1).subs{1};
                indc = L(1).subs{2};
            end
            
            % Convert ':' to actual indices
            if strcmp(indr,':')
                indr = 1:a.matdim(1);
            end
            if strcmp(indc,':')
                indc = 1:a.matdim(2);
            end
            
            % correction to allow index exceeding, expand matrix dimension
            if any(indr>a.matdim(1))
                a.matdim(1) = max(indr);
            end
            if any(indc>a.matdim(2))
                a.matdim(2) = max(indc);
            end
%             if strcmp(indc,':')
%                 a.matdim(2) = 1;
%             elseif any(indc>a.matdim(2))
%                 a.matdim(2) = max(indc);
%             end
            
            % expand coefficient matrix and reassign values 
            tmp = a.C;
            [i,j,v] = find(tmp);
            a.C = sparse(i,j,v,a.matdim(1)*(length(a.dvarname)+1),a.matdim(2)*size(a.degmat,1),nnz(a.C)+nnz(RHS.C));
            
            % To adjust these elements, we need to adjust the coefficient
            % matrix a.C. For this, identify corresponding locations of 
            % coefficients that need to be changed
            [idxI,idxJ, szR, szC] = getCindices(a,indr,indc);
            if any(size(RHS.C)~=[szR,szC])
                error('Subassignment error: dimensions of objects on the left and the right of the equation do not match');
            end
            % Then, adjust the coefficients as specified
            a.C(idxI,idxJ) = RHS.C;
        else % for now this is unused
            disp("subassgn of the type 'dpvar(a,b).prop = RHS;' is currently unsupported");
        end
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



function [indr,indc,szR,szC] = getCindices(Dp,indr,indc)
% this converts row/col indices to corresponding coeff indices in dpvar and
% returns number of rows/cols being modified
Zdlen = length(Dp.dvarname) + 1;
Zplen = size(Dp.degmat,1);

if strcmp(indr,':')
    indr = 1:size(Dp.C,1);
    szR = size(Dp.C,1);
elseif length(indr)>=2
    mdim = length(indr);
    indr = kron(vec((indr-1)*Zdlen),ones(Zdlen,1)) + repmat((1:Zdlen)',mdim,1);
    szR = length(indr);
else
    indr = (indr(1)-1)*Zdlen+1:(indr(1))*Zdlen;
    szR = length(indr);
    %szR = 1;
end

if strcmp(indc,':')
    indc = 1:size(Dp.C,2);
    szC = size(Dp.C,2);
elseif length(indc)>=2
    mdim = length(indc);
    indc = kron(vec((indc-1)*Zplen),ones(Zplen,1)) + repmat((1:Zplen)',mdim,1);
    szC = length(indc);
else
    indc = (indc(1)-1)*(Zplen)+1:(indc(1))*(Zplen);
    szC = length(indc);
    %szC = 1;
end

if any(indr<=0) || any(indc<=0)
    error('Indices must be positive integers');
end
end

function a = col_row_remove(a,L)
% Remove columns and rows of a specified by L.
% Note that L must refer to full rows or columns of a, it cannot refer to
% just a single element a(i,j) of a.

if length(L(1).subs)==1 && strcmp(L(1).subs{1},':')
    a = dpvar([]);
    return
elseif length(L(1).subs)==1 && a.matdim(1)==1
    indr = 1;
    indc = L(1).subs{1};
elseif length(L(1).subs)==1
    indr = L(1).subs{1};
    indc = ':';
else
    indr = L(1).subs{1};
    indc = L(1).subs{2};
end

if strcmp(indr,':')
    if strcmp(indc,':')
        a = dpvar([]);
    else
        cfull = 1:a.matdim(2);
        cretain = setdiff(cfull,indc);
        Ctemp = a.C;
        nz = size(a.degmat,1);
        cCretain = cell2mat(cellfun(@(x) (x-1)*nz+1:x*nz,num2cell(cretain),'uni',0));
        a.C = Ctemp(:,cCretain);
        a.matdim = [a.matdim(1),length(cretain)];
        a = combine(a);
    end
elseif strcmp(indc,':')
    rfull = 1:a.matdim(1);
    rretain = setdiff(rfull,indr);
    Ctemp = a.C;
    nd = length(a.dvarname)+1;
    rCretain = cell2mat(cellfun(@(x) (x-1)*nd+1:x*nd,num2cell(rretain),'uni',0));
    a.C = Ctemp(rCretain,:);
    a.matdim = [length(rretain),a.matdim(2)];
    a = combine(a);
else
    error('Null assignment must be specified for entire row or column')
end

end