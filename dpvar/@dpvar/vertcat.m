function G = vertcat(varargin)
% G = vertcat(E,F) performs the vertical concatenation operation [E; F]
% 
% INPUTS:
% E,F: dpvar, polynomial, or double objects (at least one is a dpvar object)
%       - Must have the same number of columns
%       - polynomials will be converted to dpvars using poly2dpvar
%       - doubles will be converted to constant dpvars
% 
% OUTPUTS:
% G: dpvar object
%       - G = [E; F]
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
% Initial coding DJ, MP, SS - 09/27/2021

if nargin==1
    % In case of just one input, no concatenation is needed
    G = varargin{1};
else
    % Otherwise, start with concatenation of first two inputs
    E = varargin{1}; F = varargin{2};
    if all(size(E)==0)
        G = F;
        return
    elseif all(size(F)==0)
        G = E;
        return
    elseif ~(size(E,2)==size(F,2))
        error('Objects being concatenated vertically have different numbers of columns');
    end
    
    % If input is not a dpvar, convert to dpvar object
    if isa(E,'double')
        E = dpvar(E,zeros(1,0),{},{},size(E));
    elseif isa(F,'double')
        F = dpvar(F,zeros(1,0),{},{},size(F));
    elseif isa(E,'polynomial')
        E = poly2dpvar(E);
    elseif isa(F,'polynomial')
        F = poly2dpvar(F);
    elseif ~isa(E,'dpvar') || ~isa(F,'dpvar')
        error('Concatenation is only supported for objects of type double, polynomial, or dpvar');
    end
    
    % To concatenate, adjust objects so that they have the same degmats,
    % dvarnames, and independent varnames
    [E,F] = common_basis(E,F);
    
    % All that remains then is concatenating the coefficient matrices, and 
    % adjusting the row dimension
    C = [E.C; F.C];
    mdim = [E.matdim(1)+F.matdim(1), E.matdim(2)];
    G = dpvar(C,E.degmat,E.varname,E.dvarname,mdim);
    
    % In case of additional inputs, simply repeat the process
    if nargin>2
        G = vertcat(G,varargin{3:end});
    end
end
end