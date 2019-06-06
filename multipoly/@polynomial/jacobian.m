function J = jacobian(f,x)
% function J = jacobian(F,X)
%
% DESCRIPTION
%   Compute the Jacobian matrix. The (i,j)-th entry of J is dF(i)/dX(j).
%
% INPUTS
%   F: Polynomial to differentiate (N-by-1 polynomial)
%   X: Variable for differentiation (V-by-1 vector of pvars
%        or cell array of strings)
%
% OUTPUTS
%   J: Jacobian of F with respect to X (N-by-V polynomial)
%
% SYNTAX
%   J = jacobian(F);
%     Computes the Jacobian of F with respect to F.varname
%   J = jacobian(F,X);
%     Computes the Jacobian of F with respect to X
%
% EXAMPLE
%   pvar x y z;
%   f = [x^3+5*y*z; 2*x*z; 3*x+4*y+6*z];
%   J = jacobian(f,[x;y;z])
%
% See also: diff, int

%   PJS 4/30/06  Initial coding for polynomial objects
%   PJS 5/21/09  Removed calls to subsasgn and diff


if nargin == 1
    x = f.varname;
end

% Convert x to a cell array of strings
if ispvar(x)
    x = char(x);
elseif ischar(x)
    x = cellstr(x);
elseif ~iscellstr(x)
    error('X must be an array of polynomial variables or strings');
end
x = x(:);
Nx = size(x,1);

% Get polynomial info about F
fcoef = f.coefficient;
fdeg = f.degmat;
fvar = f.varname;
szf = f.matdim;
nrf = szf(1);
ncf = szf(2);
Nf = nrf*ncf;

%OLD VERSION
%--------------------
% % Compute Jacobian  %Commented out by PJG 15/10/2013
% Jcoef = sparse([]);
% Jdeg = sparse([]);
% for i1=1:Nx
%     xi = x(i1);
%     
%     % Find variable we are differentiating with respect to.
%     varnumb = find( strcmp(fvar,xi) );
%     
%     % Differentiate
%     if isempty(varnumb)
%         Jcoef = [Jcoef sparse(size(Jcoef,1),Nf)];
%     else
%         tmpcoef = fcoef.*repmat(fdeg(:,varnumb),[1 Nf]);
%         Jcoef = blkdiag(Jcoef,tmpcoef);
%         
%         tmpdeg = fdeg;
%         tmpdeg(:,varnumb) = max( fdeg(:,varnumb)-1 , 0);
%         Jdeg = [Jdeg;tmpdeg];
%     end
%     
% end
%--------------------

%Compute Jacobian : PJG version, 15/10/2013
%--------------------
Jcoef = sparse(0,Nx*Nf);
Jdeg  = sparse([]);

%find the locations of any x in fvar
[tf,loc] = ismember(x,fvar);
xidx     = find(tf);

for i1 = xidx(:)'
    
    varnumb = loc(i1);
    tmprow  = sparse(size(fdeg,1),Nf*Nx);
    tmpcoef = fcoef.*repmat(fdeg(:,varnumb),[1 Nf]);
    tmprow(:,((i1-1)*Nf+1):(i1*Nf)) = tmpcoef;
    Jcoef = [Jcoef;tmprow];
    
    tmpdeg = fdeg;
    tmpdeg(:,varnumb) = max( fdeg(:,varnumb)-1 , 0);
    Jdeg = [Jdeg;tmpdeg];
end
%--------------------



chkval = 0; % skip validity check
J = polynomial(Jcoef,Jdeg,fvar,[Nf Nx],chkval);
J = combine(J);
