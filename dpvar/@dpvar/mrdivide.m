function G = mrdivide(E,F)
% G = mrdivide(E,F) performs the right matrix-division operation E/f 
% 
% INPUTS:
%   E: dpvar lass object
%   F: double with same number of columns as E
% 
% OUTPUTS:
% G: a dpvar object, G = E/F, so that G*F = E (if such a G exists)
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% S. Shivakumar at sshivak8@asu.edu, or Declan Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2023  M. Peet, S. Shivakumar, D. Jagt
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
% Initial coding DJ  - 04/04/2023

% Check that the denominator is of type 'double'.
if isa(F,'polynomial') || isa(F,'dpvar')
    try 
        F = double(F);
    catch
        error('Division of ''dpvar'' is only supported for denominators of type ''double''.')
    end
elseif ~isa(F,'double')
    error('Division of ''dpvar'' is only supported for denominators of type ''double''.')
end

if all(size(F)==[1,1])
    % Perform division by a scalar.
    if F==0
        error('Division of ''dpvar'' by 0 is not supported.')
    else
        % Divide all coefficients of G by f.
        G = E;
        G.C = G.C/F;
    end
elseif size(E,2)~=size(F,2)
    error('Column dimensions of numerator and denominator must match.')
else
    % Perform division by a matrix.
    [nr,nc] = size(E);
    nz = size(E.degmat);
    G = zeros(nr,nc);

    % Decompose E into sum of Ejj, with each jj corresponding to one row of
    % E.degmat.
    for jj=1:nz
        Cjj = E.C(:,jj:nz:end)/F;
        Gjj = dpvar(Cjj,E.degmat(jj,:),E.varname,E.dvarname,[nr,nc]);
        G = G+Gjj;
    end

% % % % % % Alternative implementation (not tested and likely slower...)
%     nd = length(E.dvarname)+1;
%     Cvec = sparse(nz*nr*nd,nc);
%     for jj=1:nz
%         Cvec((jj-1)*nr*nd+1:jj*nr*nd,:) = E.C(:,(jj-1)*nc+1:jj*nc);
%     end
%     Cvec_new = Cvex/F;
%     Cnew = sparse(nr*nd,nc*nz);
%     for jj=1:nz
%         Cnew(:,(jj-1)*nc+1:jj*nc) = Cvec_new((jj-1)*nr*nd+1:jj*nr*nd,:);
%     end
%     G = dpvar(Cnew,E.degmat,E.varname,E.dvarname,[nr,nc]);
end

end