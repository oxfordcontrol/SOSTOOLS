function varargout = subsref(a,L)
% function B=subsref(A,L);
%
% DESCRIPTION
%   Subsreference for polynomial objects.
%
% INPUTS
%   A: polynomial
%   L: a structure array with the fields:
%    type -- string containing '()'  or '.' specifying the
%              subscript type.
%    subs -- Cell array or string containing the actual subscripts.
%
% OUTPUTS
%   B: object after referencing
%
% SYNTAX
%   B=subsref(A,L)

% 6/7/2002: PJS  Initial Coding
% 06/22/2022: DJ, Add option for logical row and column indexing
% 11/25/2026: DJ, Linear indexing gives the output shape MATLAB would, and
%   accepts matrix-valued indices;

switch L(1).type
    case '.'
        b = get(a,L(1).subs);
    case '()'
        sza = a.matdim;
        nra = sza(1);
        nca = sza(2);
        acoef = a.coefficient;
        b = a;
        
        if length(L(1).subs)==1
            % Check Indices
            subsidx = L(1).subs{1};
                
            if strcmp(subsidx,':')
                subsidx = (1:nra*nca)';
            end
            if islogical(subsidx)
                subsidx = find(subsidx);
            end
            
            if ~isempty(subsidx) && ...
                    (max(subsidx(:))>nra*nca || min(subsidx(:))<1)          % DJ, 11/25/2026
                error('Index exceeds matrix dimensions.');
            end
            
            % Do subsref
            b.coefficient = acoef(:,subsidx);
            if nra*nca>1 && min(nra,nca)==1 && isvector(subsidx) && ~strcmp(L(1).subs{1},':')    % DJ, 11/25/2026
                % A vector indexed by a vector keeps its orientation
                b.matdim = [length(subsidx) 1];
                if nra==1,  b.matdim = [1 length(subsidx)];  end
            else
                % Otherwise the output takes the shape of the index
                b.matdim = size(subsidx);
            end
    
        elseif length(L(1).subs)==2
            ridx = L(1).subs{1};
            if strcmp(ridx,':')
                ridx = 1:nra;
            end
            if islogical(ridx)  % DJ, 06/22/2022
                ridx = find(ridx);
            end
            cidx = L(1).subs{2};
            if strcmp(cidx,':')
                cidx = 1:nca;
            end
            if islogical(cidx)
                cidx = find(cidx);
            end
            
            % Check Indices
            if ~isempty(ridx) && (max(ridx)>nra || min(ridx)<1)
                error('Index exceeds matrix dimensions.');
            elseif ~isempty(cidx) && (max(cidx)>nca || min(cidx)<1)
                error('Index exceeds matrix dimensions.');
            end
            
            % Do subsref
            bcoef = acoef;
            idx = [];
            for i1 = 1:length(cidx)
                idx = [idx ridx+(cidx(i1)-1)*nra];
            end
            b.coefficient = bcoef(:,idx);
            b.matdim = [length(ridx) length(cidx)];
            
        else
            error('Invalid subsref for polynomials');
        end
        
%         % Remove any terms with zero coeffs and then combine down
%         % (This is faster than directly calling combine without
%         %  pruning out terms with zero coefs)
%         bcoef = b.coefficient;
%         [ridx,cidx] = find( bcoef );
%         ridx = unique(ridx);
%         b.coefficient = bcoef(ridx,:);
%         b.degmat = b.degmat(ridx,:);
        
        if isempty(b.coefficient)
            b = polynomial(zeros(b.matdim));
        else
            b = combine(b);
        end
        
    case '{}'
        error('{}- like subsreference is not supported for polynomial objects.');
end

if length(L)>1
    b = subsref(b,L(2:end));
end
varargout{1}=b;


