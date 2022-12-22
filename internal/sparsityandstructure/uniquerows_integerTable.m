function [newC,uniqueZ] = uniquerows_integerTable(C,Z,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [newC,uniqueZ] = uniquerows_integerTable(C,Z,opts) removes
% repeated rows from integer array Z, and combines the associated
% columns in C such that
%   newC * uniqueZ = C * Z;
%
% INPUTS:
% - C:      A m x n array of coefficients, sparse or not.
%           If opts=true or opts='transpose', will by of dim n x m.
% - Z:      A n x p array of nonnegative integers.
% - opts:   Optional string 'transpose', or logical true/false. If 
%           specified and true, or set to 'transpose', output will satisfy
%               newC' * uniqueZ = C' * Z;
% 
% OUTPUTS:
% - newC:       A m x k array of coefficients, sparse if input C is sparse.
%               If opts=true or opts='transpose', will be of dim k x m.
% - uniqueZ:    A k x p array of nonnegative integers, with each row equal
%               to some row in the input Z, but with no two rows the same.
%
% NOTES: 
% If no coefficients C are specified, the function assumes C = eye(nZ),
% such that   newC * uniqueZ = Z   or   newC' * uniqueZ = Z;
% If also nargout==1, only the unique rows uniqueZ will be returned, no
% coefficient matrix newC.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This file is part of SOSTOOLS - Sum of Squares Toolbox ver 4.00.
%
% Copyright (C)2002, 2004, 2013, 2016, 2018, 2021  
%                                      A. Papachristodoulou (1), J. Anderson (1),
%                                      G. Valmorbida (2), S. Prajna (3), 
%                                      P. Seiler (4), P. A. Parrilo (5),
%                                      M. Peet (6), D. Jagt (6)
% (1) Department of Engineering Science, University of Oxford, Oxford, U.K.
% (2) Laboratoire de Signaux et Systmes, CentraleSupelec, Gif sur Yvette,
%     91192, France
% (3) Control and Dynamical Systems - California Institute of Technology,
%     Pasadena, CA 91125, USA.
% (4) Aerospace and Engineering Mechanics Department, University of
%     Minnesota, Minneapolis, MN 55455-0153, USA.
% (5) Laboratory for Information and Decision Systems, M.I.T.,
%     Massachusetts, MA 02139-4307
% (6) Cybernetic Systems and Controls Laboratory, Arizona State University,
%     Tempe, AZ 85287-6106, USA.
%
% Send bug reports and feedback to: sostools@cds.caltech.edu
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
% along with this program. If not, see <http://www.gnu.org/licenses/>.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change log and developer notes
% 12/17/2022 -- DJ: Copy from @dpvar/private/DPVuniqueterm
%

% Check the inputs
use_eye = false;        % use_eye = true        -->  C = I;
use_transpose = false;  % use_transpose = true  -->  newC'*uniqueZ = C'*Z
if nargin==1
    % Assume C=I, and first input corresponds to Z.
    use_eye = true;
    Z = C;
elseif nargin==2 && (ischar(Z) || islogical(Z))
    % Assume C=I, first input corresponds to Z, and second to opts.
    use_eye = true;
    opts = Z;
    Z = C;
    % Check the opts are properly specified.
    if islogical(opts)
        use_transpose = opts;
    elseif ischar(opts) && strcmpi(opts,'transpose')
        use_transpose = true;
    else
        error('Second argument is not appropriate. Call "help uniquerows_integerTable" for more information.')
    end
elseif nargin==3
    % Check the opts are properly specified.
    if islogical(opts)
        use_transpose = opts;
    elseif ischar(opts) && strcmpi(opts,'transpose')
        use_transpose = true;
    else
        error('Third argument must be one of ''transpose'' (char) or true/false (logical).')
    end
elseif nargin>3
    error('Too many inputs.')
end

% Find repeated monomials
% Sort done by total total monomial degree first and then lexicographic
%   sortidx: index such that degsort = degmat(sortidx)
%   repeatidx: nonzero indices indicate repeats
%   rvec:  map rows of degsort to numbers 1,..,nut
M = [sum(Z,2) Z];
[degsort,sortidx] = sortrows_integerTable(M);
degsort = degsort(:,2:end);
repeatidx = ~any(degsort(2:end,:) ~= degsort(1:end-1,:),2);
rvec = cumsum([1; ~repeatidx]);

% Create unique set of monomials
uniqueidx = sortidx;
uniqueidx( repeatidx ) = [];
uniqueZ = Z(uniqueidx,:);
if nargout==1 && use_eye
    newC = uniqueZ;
    return
end

% Sum repeated coefficients in C
%  (Use sparse multiply, C*summat, to sum the columns)
nt = size(Z,1);
nut = size(uniqueZ,1);
if ~use_transpose
    % C * Z = C * (summat*uniqueZ) = newC * uniqueZ;
    summat = sparse(sortidx,rvec,1,nt,nut,nt);
    if use_eye
        % C = I;
        newC = summat;
    else
        newC = C*summat;
    end
else
    % C' * Z = C' * (summat'*uniqueZ) = (summat*C)' * uniqueZ = newC' * uniqueZ;
    summatT = sparse(rvec,sortidx,1,nut,nt,nt);
    if use_eye
        % C = I;
        newC = summatT;
    else
        newC = summatT*C;
    end
end

end