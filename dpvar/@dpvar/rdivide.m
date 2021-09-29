function E = rdivide(E,f)
% G = rdivide(E,f) performs the right division operation E./f 
% 
% INPUTS:
%   E: dpvar or polynomials class object
%   f: double
% 
% OUTPUTS:
% G: a dpvar object, G = E./f
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
% Initial coding DJ, MP, SS  - 09/27/2021

% Check input types
if ~isa(f,'double')
    try
        f = double(f);
    catch
        error('For division involving dpvars, denominator must be a non-zero double')
    end
end
if any(f==0,'all') || (any(size(f)~=size(E)) && (numel(f)~=1))
    error("For division involving dpvars, denominator must be non-zero scalar or double of same size as numerator");
end

if numel(f)==1
    % Scalar division
    E.C = E.C./f;
else
    % Divide coefficients according to the matrix element they're associated with
    nd = length(E.dvarname)+1;
    nz = size(E.degmat,1);

    F = kron(f,ones(nd,nz));
    E.C = E.C./F;
end

end