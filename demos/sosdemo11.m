% SOSDEMO11 --- Robust Stability
%
% We have a linear system of the form
% \dot x(t) = A(p)x(t) 
% where we would like to prove stability for all p \in G:=\{p: g(p)\ge 0\}
% To this end, we
% Construct P(p) > \epsilon I           for all p \in G
% such that A(p)^T P(p)+P(p)A(p)<0      for all p \in G
%
clear; echo on;
pvar p1 p2
vartable = [p1 p2];

eps = .01;

% =============================================
% This is the problem data
n = 3;      % Size of matrix A
deg = 3;    % Monomial degree for solving the program

% =============================================
% Set up the uncertain matrix, \dot{x}=A(p)*x
A = -50*eye(n) + .05*tril(ones(n))*p1 - .05*triu(ones(n))*p2;

% =============================================
% Set up the domain constraint, g(p)>=0
g = (1-p1^2-p2^2); 

% =============================================
% Initialize the sum of squares program
prog = sosprogram(vartable);

% =============================================
% Define the positive matrix P
Zsym = monomials(vartable,0:2*deg);
[prog,P] = sospolymatrixvar(prog,Zsym,[n n],'symmetric');
[prog] = sosmatrixineq(prog,P-eps*eye(n));

% =============================================
% Define the SOS multiplier
[prog,s] = sospolymatrixvar(prog,Zsym,[n n],'symmetric');
[prog] = sosmatrixineq(prog,s);

% =============================================
% Define the derivative matrix A'*P+P*A
Dif = A'*P + P*A;

% =============================================
% Require a negative derivative in the domain
[prog] = sosmatrixineq(prog,-s*g-Dif);

% =============================================
% Solve the problem
solver_opt.simplify = 'on';

solver_opt.solver = 'sedumi';
prog = sossolve(prog,solver_opt);

% If program is feasible, the system \dot{x}=A(p)*x is stable for any
% p such that p'*p<=1 