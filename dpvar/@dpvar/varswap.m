function E = varswap(E, pvars1, pvars2)
% E = varswap(E, pvars1, pvars2) swaps the independent variables pvars1 
% with pvars2 in the dpvar E
% 
% INPUTS:
% E: dpvar
% pvars1: polynomial or cell of size nx1
%       - if polynomial, each element must be a single polynomial variable
%       - if cell, each element must be a char specifying a variable name
% pvars2: object of same type and size as pvars1
% 
% OUTPUTS:
% E: dpvar
%       - E.varname will be adjusted so that, for each i=1:n, pvars1{i} and 
%           pvars2{i} are exchanged 
%
% EXAMPLE:
% Let D(d1,d2,x1,x2,x3,x4) be a dpvar with decision vars d1, d2, and
% independent vars x1, x2, x3, x4. Then:
%   - E = varswap(D,{'x1';'x3'},{'x2';'x4'});
%   - E = varswap(D,{'x3';'x1'},{'x4';'x2'});
%   - E = varswap(D,[x1;x3],[x4;x2]); with [x1;x3] and [x4;x2] polynomial
% will all produce the same dpvar E, so that
%   E(d1,d2,x1,x2,x3,x4) = D(d1,d2,x2,x4,x1,x3)
%
% NOTES:
% - To avoid issues, make sure that:
%   > pvars1 and pvars2 do not contain any of the same variables
%   > pvars1 and pvars2 both contain no repeats
%   This is to avoid swapping one variable twice (x with y and x with z),
%   for which there may be multiple options
%
% - This function can NOT be used to swap decision variables, only
%   independent variables
%
% - If an an element pvars1{i} does not appear as variable in E, but the 
%   associated element pvars2{i} does appear in E, the name of this 
%   variable will be changed to that of pvars1{i}.
%   Conversely, if an an element pvars2{i} does not appear in E but 
%   pvars1{i} does, the name of this variable will be changed to that of
%   pvars2{i}. See also the function "subs" for variable substitution.
%   If neither variable appears, these elements will be ignored.
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
% Initial coding DJ, MP, SS - 07/12/2021

% In case input pvars1 is polynomial, convert it to a cell of varnames
if isa(pvars1,'polynomial')
    tmp = cell(1,length(pvars1));
    for i=1:length(pvars1)
        tmp{i} = pvars1(i).varname{1};
    end
    pvars1 = tmp;
elseif ~isa(pvars1,'cell')
    pvars1 = {pvars1};
end
% In case input pvars2 is polynomial, convert it to a cell of varnames
if isa(pvars2,'polynomial')
    tmp = cell(1,length(pvars2));
    for i=1:length(pvars2)
        tmp{i} = pvars2(i).varname{1};
    end
    pvars2 = tmp;
elseif ~isa(pvars2,'cell')
    pvars2 = {pvars2};
end

% Make sure the sizes of the varnames match
if length(pvars1)~=length(pvars2)
    error('Length of old variable names and new variable names must be equal');
end

% Let vlist be the old varnames, and tmp the adjusted varnames
vlist = E.varname;
tmp = vlist;

% Iteratively swap the independent variable names as specified
for i=1:length(pvars1)
    oldidx = find(strcmp(vlist,pvars1{i}));     newidx = find(strcmp(vlist,pvars2{i}));
    if ~isempty(oldidx)
        tmp{oldidx} = pvars2{i};
    end
    if ~isempty(newidx)
        tmp{newidx} = pvars1{i};
    end
end

% Assign the new varnames to E
E.varname = tmp;
end