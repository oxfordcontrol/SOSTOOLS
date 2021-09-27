% SOSDEMO4s --- Matrix Copositivity
% Section 4.4 of SOSTOOLS User's Manual

clear; echo on;
syms x1 x2 x3 x4 x5;
vartable = [x1; x2; x3; x4; x5];

% The matrix under consideration
J = [1 -1  1  1 -1;
    -1  1 -1  1  1;
     1 -1  1 -1  1;
     1  1 -1  1 -1;
    -1  1  1 -1  1];

% =============================================
% First, initialize the sum of squares program
prog = sosprogram(vartable);     % No decision variables.

% =============================================
% Next, define SOSP constraints

% Constraint : r(x)*J(x) - p(x) = 0
J = [x1^2 x2^2 x3^2 x4^2 x5^2]*J*[x1^2; x2^2; x3^2; x4^2; x5^2];
r = x1^2 + x2^2 + x3^2 + x4^2 + x5^2;

prog = sosineq(prog,r*J);

% =============================================
% And call solver
solver_opt.solver = 'sedumi';
prog = sossolve(prog,solver_opt);

% =============================================
% If program is feasible, the matrix J is copositive.
echo off