function varargout = size(A,dim)
% varargout = size(A,dim) returns the size of a dpvar along dimension dim
% 
% INPUTS:
% A: dpvar object
% dim: scalar value 1 or 2 (optional input)
%       - if 1, number of rows of A is returned as output
%       - if 2, number of columns of A is returned as output
%       - if not specified, full matdim of A will be returned
% 
% OUTPUTS:
% out: scalar or 1x2 array
%       - if one output is called and dim is specified, the size of the
%           dpvar matrix along the specified dimension is returned
%       - if one output is called and no dim is specified, a 1x2 array
%           specifying [row,col] dimension of the dpvar is returned
%       - if two outputs are called, the first output will be assigned the
%           number of rows, and the second the number of columns
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
% Initial coding, DJ, MP, SS - 06/10/2021
  
if nargin ==1 % if no dim is specified, return both row and column numbers 
    out = A.matdim;
    if nargout==0 || nargout==1
        varargout{1}=out;
    elseif nargout==2   % allow row and column numbers to be output separately
        varargout{1} = out(1);
        varargout{2} = out(2);
    end 
elseif nargin == 2
    if dim==1 || dim==2
        out = A.matdim(dim);
    else
        error('Dimension must be 1 or 2');
    end
    varargout{1} = out;  
end
end