% SOSDEMO12s --- Local Stability
%
% Testing local stability of coupled Van der Pol equations, \dot{x}=f(x),
% on a circle x*x'<=r. In particular, x=[y,z] and:
% % dot(yi) = -2*zi                                         % i=1 to N
% % dot(zj) = 0.8*yj + 10*(1.2^2*yj^2-0.21)*zj+ep_jz(j+1)yj % j=1 to N-1
% % dot(zN) = 0.8*yN + 10*(1.2^2*yN^2-0.21)*zN
%
% See also
% Tacchi et al., Approximating Regions of Attraction of a Sparse Polynomial
% Differential System, 2020.
%
clear; echo on;

eps = 0.001;

% =============================================
% This is the problem data
n = 3;                  % Number of coupled ODEs
deg = 3;                % Monomial degree for solving the program
r = 0.5;                % Radius of circle in which to test stability
ep_c = randn(1)*0.5-1;  % Function parameter


try 
% =============================================
% Set up independent variables
svartable = sym(zeros(1,2*n));
for k=1:n
    var_y = ['y',int2str(k)];
    var_z = ['z',int2str(k)];
    svartable(2*k-1) = sym(var_y,'real');
    svartable(2*k) = sym(var_z,'real');
end

% =============================================
% Set up the nonlinear function, \dot{x}=f(x)
f_new = sym(zeros(2*n,1));
for m=1:2:2*n
    f_new(m) = -2*svartable(m+1);    % xvartable = [y1, z1, y2, z2,..., yK,zK]
end
for m=2:2:2*n-1
    f_new(m) = 0.8*svartable(m-1) + 10*(1.2^2*svartable(m-1)^2-0.21)*svartable(m)+ep_c*svartable(m+2)*svartable(m-1);
end
f_new(2*n) = 0.8*svartable(end-1) + 10*(1.2^2*svartable(end-1)^2-0.21)*svartable(end);

% =============================================
% Set up the domain constraint, g(x)>=0
g = r - (svartable*svartable');

% =============================================
% Initialize the sum of squares program
prog = sosprogram(svartable);

% =============================================
% Define the quadratic Lyapunov function
Z = monomials(svartable,1:deg);
[prog,V] = sossosvar(prog,Z);
V = V + eps*(svartable*svartable');

% =============================================
% Define the SOS multiplier
[prog,s] = sossosvar(prog,Z);

% =============================================
% Define the derivative of the Lyapunov Function
Vd = jacobian(V,svartable)*f_new;

% =============================================
% Require a negative derivative in the domain
[prog] = sosineq(prog,-Vd-s*g);

% =============================================
% Solve the problem
solver_opt.solver = 'sedumi';
prog = sossolve(prog,solver_opt);

% If program is feasible, the system \dot{x}=f(x) is stable in the
% circle x^x' <= r.

catch
    
disp('Problem may be too big to handle with syms. Please use a smaller problem, or try the ''dpvar'' implementation')
    
end

echo off;
