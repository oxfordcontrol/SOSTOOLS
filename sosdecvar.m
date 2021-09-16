function sos = sosdecvar(sos,decvartable)
% SOSDECVAR --- Declare new decision variables in an SOS 
%       program.
%
% SOSP = sosdecvar(SOSP,DECVARS)
%
% SOSP is the sum of squares program.
% DECVARS is a vector whose entries are the new decision 
% variables. DECVARS is either a symbolic or a polynomial object.
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
% 20/03/02 - SP
% 07/28/21 - MP Conversion to dpvar

% Update tables and indices
if isfield(sos,'symvartable')

	decvartabletemp = sym2chartable(decvartable);  
	if length(sos.decvartable) ~= 2
        sos.decvartable = [decvartabletemp(1:end-1),',',sos.decvartable(2:end)];
        
    else
        sos.decvartable = [decvartabletemp(1:end)];  % Special case: previously empty table
	end;
	if size(decvartable,2) > 1
        decvartable = decvartable.';
	end;
	sos.symdecvartable = [decvartable; sos.symdecvartable];
    offset = length(decvartable);
elseif isa(decvartable,'dpvar')
    % this puts all decision variables in the DPvar into the prog decision
    % variable table.
    sos.decvartable = [decvartable.dvarname; sos.decvartable];
    offset = length(decvartable.dvarname);
else
    % This converts all variables in the poly to dpvars and 
    % puts all these variables into the prog decision
    % variable table.
    
    sos.decvartable = [decvartable.varname; sos.decvartable];
    % The decision was made not to promote pvars to dpvars in order to
    % avoid conflicts with previous pvar structures containing the decision
    % variable names
%    dpvar_out = poly2dpvar(decvartable, decvartable.varname);
%     for i=1:length(decvartable)
%     assignin('base',decvartable(i).varname{1},dpvar_out(i))
%     end
    offset = length(decvartable.varname);
end;
    
    
%offset = length(decvartable); -MP moved this into the cases to reflect
%different uses for poly/dpvar/sym
for i = 1:sos.var.num+1
    sos.var.idx{i} = sos.var.idx{i} + offset;
end;

% Update existing equations
for i = 1:sos.expr.num
    sos.expr.At{i} = [sparse(offset,size(sos.expr.At{i},2));sos.expr.At{i}]; 
end;

% Update existing objective
sos.objective = [sparse(offset,1); sos.objective];


