function [xt,ft,flg] = psolve(f,x0,opts)
% function [xt,ft,flg] = psolve(f,x0,opts)
%
% DESCRIPTION 
%   This function solves for roots of the polynomial equation
%      f(x)=0
%   FSOLVE is used to solve these nonlinear equations.  
%   Initial guesses for the trim state/input can be passed to FSOLVE.  
%
% INPUTS 
%   f: System of polynomial equations  (Nx-by-1 polynomial)
%   x0: Initial guess for solution [Optional, Default: x0=0]
%   opts: options structures may be passed to fsolve [Optional]
%     
% OUTPUTS
%   xt: Trim state (Nx-by-1 vector)
%   ft: f evaluated at (xt) (Nx-by-1 vector)
%      If psolve was successful finding a trim point then ft:=f(xt,ut)
%      will be equal to zero
%   flg: Exit flag returned by fsolve
%
% SYNTAX
%   [xt,ft,flg] = psolve(f)
%   [xt,ft,flg] = psolve(f,x0)
%   [xt,ft,flg] = psolve(f,x0,opts)
%
% EXAMPLE
%   pvar x1 x2 u;
%   f = [-2*x1+x2+x1^2-7; x1-3*x2+3];
%
%   % Find a root of the polynomial system
%   [xt,ft] = psolve(f)
%
%   % Find a root using opts structure
%   x0 = []; opts=optimset('Display','On');
%   [xt,ft] = psolve(f,x0,opts)
%
% See also fsolve, plinearize

% MMP 9/13/2021   Initial Coding

% Error Checking
if nargin==1
    x0 = [];
    opts = [];
elseif nargin ==2
    opts = [];        
end

% Set default initial starting guesses 
Nx = length(f.varname);
if isempty(x0)
    x0 = zeros(Nx,1);
end

% Set default fsolve options
if isempty(opts)
    opts = optimset('Display','Off');    
end

% Identify variable names
x=polynomial(f.varname);

% Solve for roots, z:=[x;u]
fh = @(z) double(subs(f,x,z)); 
[xt,ft,flg] = fsolve(fh,x0,opts);

