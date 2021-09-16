function polyVar = dpvar2poly(dpVar)
% This function converts an object of dpvar class to polynomial class
% object. Function can also take cell arrays as inputs where all cells are
% dpvar class objects.
% 
% INPUTS:
% dpVar: a dpvar class object or a cell array of dpvar class objects. Use
% "help dpvar" command to check the structure
% 
% OUTPUTS:
% polyVar: polynomial class equivalent of dpVar when input is a dpvar, else
% cell array of polynomial class objects when input is a cell array of
% dpVars. use "help polynomial" command to check the structure
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
% 
% Initial coding, DJ, MP, SS - 06/10/2021

% check if input object is dpvar or cell array of dpvars
if ~isa(dpVar, 'dpvar') && (isa(dpVar,'cell')&& ~isa(dpVar{1},'dpvar'))
    error('Input object must be of type dpvar or a cell array with dpvar objects');
end

if isa(dpVar, 'dpvar')
    if dpVar.chkval

        % matdim remains same
        mdim = dpVar.matdim;

        % put varnames and dvarnames together
        vnames = [dpVar.dvarname; dpVar.varname];

        % construct degmat to include monomials from both varnames and
        % dvarnames. dvarnames will only have degrees 1 and hence degmat is an identity matrix
        dmat_d = [[zeros(1,length(dpVar.dvarname));speye(length(dpVar.dvarname))],...
                   zeros(length(dpVar.dvarname)+1,length(dpVar.varname))];
        dmat_p = [zeros(size(dpVar.degmat,1),length(dpVar.dvarname)),dpVar.degmat];

        % now take all permutations of rows in degmat_d and degmat_p and
        % add them
        dmat = kron(dmat_p,ones(size(dmat_d,1),1)) + repmat(dmat_d,size(dmat_p,1),1);

        % reshape coeff matrix by splitting into block matrix stucture using mat2cell
        % [b11 b12 b13]  -> [b11' b12' b13']
        % [b21 b22 b23]     [b21' b22' b23']
        coeff = mat2cell(dpVar.C,(length(dpVar.dvarname)+1)*ones(size(dpVar,1),1),ones(size(dpVar.C,2),1));
        coeff = cell2mat(cellfun(@transpose,coeff,'un',0));
        if all(size(coeff)==0) % correction if row or col dimension is zero
            coeff = sparse(size(dpVar,1),size(dmat,1)*size(dpVar,2));
        end
        coeff = mat2cell(coeff,ones(size(dpVar,1),1),size(dmat,1)*ones(size(dpVar,2),1));
        coeff = cell2mat(coeff(:));
        if all(size(coeff)==0) % correction if row or col dimension is zero
            coeff = sparse(prod(mdim),size(dmat,1));
        end

        polyVar = polynomial(coeff', dmat, vnames, mdim);
    else
        error('Input object has inconsistent dimensions');
    end
elseif isa(dpVar,'cell')
    polyVar = cellfun(@(x) dpvar2poly(x), dpVar,'un',0);
end

end