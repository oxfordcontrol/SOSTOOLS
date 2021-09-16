function out = get(obj,prop)
% out = get(obj,prop) returns field prop from a dpvar object obj
% 
% INPUTS:
% obj: dpvar object
% prop: struct of length 1 or 2, with fields "type" and "subs"
%   - prop(1).type must be '.', with prop(1).subs one of the dpvar fieldnames
%     "matdim", "degmat", "C", "varname", "dvarname", or "chkval"
%   - prop(2).type can be '()' or '{}' extracting array or cell elements 
%     prop(2).subs of the field specified in prop(1)
% 
% OUTPUTS:
% out: cell, array or char
%   - Output corresponds to the desired (elements of) obj.(prop(1).subs)
% 
% EXAMPLE:
% Let Dp be a dpvar specifying a polynomial
% polynomial(Dp) = [d1, x; d2*y, x*y];
% Then:
%   - get(Dp,prop) = Dp.dvarname = {'d1';'d2'}
%           when prop.type = '.' and prop.subs = 'dvarname'
%   - get(Dp,prop) = Dp.dvarname{2} = 'd2'
%           when prop(1).type = '.' and prop(1).subs = 'dvarname' and
%                prop(2).type = '{}' and prop(2).subs = 2
%   - get(Dp,prop) = Dp.matdim = [2,2]
%           when prop.type = '.' and prop.subs = 'matdim'
%   - get(Dp,prop) = Dp.matdim(1) = 2
%           when prop(1).type = '.' and prop(1).subs = 'matdim' and
%                prop(2).type = '()' and prop(2).subs = 1
%
% NOTES:
% - In general, get will only be called by "subsref.m".
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
% Initial coding, DJ, MP, SS - 07/02/2021

% Extract the desired field from the dpvar object
switch prop(1).subs
    case 'matdim'
        out = obj.matdim;
        if length(prop)==2
            % If an additional prop is specified, extract only desired
            % elements from the field
            out = out(prop(2).subs{1});
        end
    case 'degmat'
        out = obj.degmat;
        if length(prop)==2
            if length(prop(2).subs)==1
                out = out(prop(2).subs{1});
            else
                out = out(prop(2).subs{1},prop(2).subs{2});
            end
        end
    case 'C'
        out = obj.C;
        if length(prop)==2
            if length(prop(2).subs)==1
                out = out(prop(2).subs{1});
            else
                out = out(prop(2).subs{1},prop(2).subs{2});
            end
        end
    case 'varname'
        out = obj.varname;
        if length(prop)==2
            % For cell arrays, allow both '()' and '{}' type indexing
            if strcmp(prop(2).type,'()')
                out = out(prop(2).subs{1});
            elseif strcmp(prop(2).type,'{}')
                out = out{prop(2).subs{1}};
            end
        end
    case 'dvarname'
        out = obj.dvarname;
        if length(prop)==2
            % For cell arrays, allow both '()' and '{}' type indexing
            if strcmp(prop(2).type,'()')
                out = out(prop(2).subs{1});
            elseif strcmp(prop(2).type,'{}')
                out = out{prop(2).subs{1}};
            end
        end
    case 'chkval'
        out = obj.chkval;
    end
end