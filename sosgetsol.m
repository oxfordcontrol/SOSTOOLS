function p = sosgetsol(sos,V,digit)
% 
% SOSGETSOL --- Get the solution from a solved SOS program 
%
% SOL = sosgetsol(SOSP,VAR,DIGIT) 
%
% SOL is the solution from a solved sum of squares program SOSP,
% obtained through substituting all the decision variables
% in VAR by the numerical values which are the solutions to
% the corresponding semidefinite program. 
%
% The third argument DIGIT (optional) will determine the 
% accuracy of SOL in terms of the number of digits. Default 
% value is 5.
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change log and developer notes
% 12/27/01 - SP
% 02/21/02 - SP -- Symbolic polynomial
% 03/01/02 - SP -- New syntax
% 03/15/02 - SP -- Fast
% 19/06/14 - GV, extends getsol for matrix variable under pvar
% 06/09/13 - MP -- faster implementation for matrix-valued polynomials
% 06/25/20 - Sachin -- fixed bug related constant polynomials
% 04/19/22 - DJ -- Update to allow cellstr inputs
% 07/02/24 - DJ -- Add check if program has been solved.
% 01/26/25 - DJ -- Bugfix dpvar to polynomial conversion, avoid duplicate
%                   varnames.
% 05/28/26 - DJ -- Implement distinct substitution for 'dpvar' case.

if nargin == 2
    digit = 5;   % Default
end

% Check if the program has even been solved yet.
if ~isfield(sos,'solinfo') || ~isfield(sos.solinfo,'RRx') || isempty(sos.solinfo.RRx)
    error('No solution seems to have been produced: either "sossolve" has not been called, or the program was found to be infeasible.')
end

if isfield(sos,'symvartable')
    
    p = mysymsubs(V,sos.symdecvartable,sos.solinfo.RRx(1:length(sos.symdecvartable)),digit);

else
    if ischar(V)   % DJ - 04/19/22
        V = polynomial({V});
    elseif iscellstr(V)
        V = combine(polynomial(V));
    elseif isa(V,'cell')
        p = cell(size(V));
        for k=1:numel(p)
            p{k} = sosgetsol(sos,V{k},digit);
        end
        return
    elseif isa(V,'double')
        p = V;
        return
    end

    % At this point, V should be either 'dpvar' or 'polynomial
    if isa(V,'dpvar')                                                       % DJ, 05/28/26
        % % V is a dpvar object
        % Decompose the polynomial variable
        dvarname = V.dvarname;      nZ = numel(dvarname);
        varname = V.varname;        degmat = V.degmat;
        Cmat = V.C;
        [m,n] = size(V);

        % Determine which of the variables in V are known
        [~,idxdecvar1,idxdecvar2] = intersect(dvarname,sos.decvartable);
        idxvar = setdiff((1:length(dvarname))',idxdecvar1); % indices of remaining dvars
        dvar_new = dvarname(idxvar);                        % names of remaining dvars
        nZ_new = numel(idxvar);                             % number of remaining dvars
        % Get the values of known variables
        dvals = zeros(numel(idxdecvar2),1);
        dvals(idxdecvar1) = sos.solinfo.RRx(idxdecvar2);
        dvals = [1;dvals];                      % add constant term
        idxvar = [1;idxvar+1];                  % constant term always remains
        
        % Declare new coefficients given known decision variable values
        [ridcs,cidcs,vals] = find(Cmat);
        ridcs_new = [];     cidcs_new = [];     vals_new = [];
        if isempty(dvar_new)
            % All variables are replaced
            for i=1:m
                % Get coefficient in row i of the matrix-valued function V
                is_row_i = ridcs>(i-1)*(nZ+1) & ridcs<=i*(nZ+1);
                ridcs_i = ridcs(is_row_i)-(i-1)*(nZ+1);     % index of decision variable
                vals_i = vals(is_row_i);
                cidcs_i = cidcs(is_row_i);
                % Multiply coefficients with variable values
                vals_fix = dvals(ridcs_i).*vals_i;
                % Link known values to the constant monomial
                ridcs_new = [ridcs_new; ((i-1)*(nZ_new+1)+1)*ones(size(ridcs_i))];
                vals_new = [vals_new; vals_fix(:)];
                cidcs_new = [cidcs_new; cidcs_i];
            end
        else
            % Some variables are retained
            new_idcs = zeros(numel(dvarname)+1,1);
            new_idcs(idxvar) = 1:numel(idxvar);     % new dvar index associated with each retained dvar
            for i=1:m
                % Get coefficient in row i of the matrix-valued function V
                is_row_i = ridcs>(i-1)*(nZ+1) & ridcs<=i*(nZ+1);
                ridcs_i = ridcs(is_row_i)-(i-1)*(nZ+1);     % index of decision variable
                cidcs_i = cidcs(is_row_i);
                vals_i = vals(is_row_i);
                % Extract the coefficients that are multiplied with known
                % values
                isvar = ismember(ridcs_i,idxvar);
                ridcs_fix = ridcs_i(~isvar);
                cidcs_fix = cidcs_i(~isvar);
                vals_fix = vals_i(~isvar);
                % Multiply coefficients with variable values
                vals_fix = dvals(ridcs_fix).*vals_fix;
                % Link known values to the constant monomial
                ridcs_new = [ridcs_new; ((i-1)*(nZ_new+1)+1)*ones(numel(vals_fix),1)];
                cidcs_new = [cidcs_new; cidcs_fix(:)];
                vals_new = [vals_new; vals_fix(:)];
                % Add remaining coefficients back into the list
                ridcs_var = new_idcs(ridcs_i(isvar));       % account for removed dvars
                ridcs_var = ridcs_var + (i-1)*(nZ_new+1);   % account for row number in matrix
                ridcs_new = [ridcs_new; ridcs_var];
                cidcs_new = [cidcs_new; cidcs_i(isvar)];
                vals_new = [vals_new; vals_i(isvar)];
            end
        end
        % Build the new coefficient matrix
        Cnew = sparse(ridcs_new,cidcs_new,vals_new,m*(nZ_new+1),size(Cmat,2));
        % Build a dpvar object with the known values set
        p = dpvar(Cnew, degmat, varname, dvar_new, [m,n]);
        % If possible, convert to polynomial
        if isempty(dvar_new)
            p = combine(dpvar2poly(p));
        else
            p = combine(p);
        end
    else
        % % V is polynomial
        [~,idxdecvar1,idxdecvar2] = intersect(V.varname,sos.decvartable);
        idxvar = setdiff(1:length(V.varname),idxdecvar1);
        coeffs = [V.degmat(:,idxdecvar1) 1-sum(V.degmat(:,idxdecvar1),2)]*[sos.solinfo.RRx(idxdecvar2);1];
                                                    %MMP 6/9/2013 updated to
                                                    %allow for matrix-valued
                                                    %polynomials and to fix
                                                    %problem with terms with no
                                                    %decision variables
        coeffs = V.coefficient.*repmat(coeffs,1,size(V.coefficient,2));         % 01/31/02
        varname = V.varname(idxvar);
        if isempty(idxvar)
            degmat = [];
        else
            degmat = V.degmat(:,idxvar);
        end
        if isempty(degmat)
            coeffs=sum(coeffs,1); % modified by sachin - 6/25/2020 original version sum(coeffs)
            p=polynomial(coeffs);
            p=reshape(p,size(V));
        else
            p=polynomial(coeffs,degmat,varname,size(V));
        end
        p=combine(p);
    end
end
