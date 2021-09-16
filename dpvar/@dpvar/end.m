function out = end(a,k,n)
% out = end(a,k,n) is used when extracting elements of a dpvar object up to
% its end, e.g. a(3:end,1) or a(:,end)
% 
% INPUTS: 
%  There are no inputs. Using 'end' while indexing automatically calls this function
% 
% OUTPUTS: 
% Index corresponding to last row/col
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
% DJ, MP, SS, 07/28/2021: Copied from "@polynomial" implementation

if n==1
    sza = size(a);
    out = sza(1)*sza(2);
elseif n==2
    out = size(a,k);
end
end