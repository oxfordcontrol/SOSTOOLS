% SOSDEMO1s --- Sum of Squares Test
% Section 4.1 of SOSTOOLS User's Manual

clear; echo on;
syms x1 x2;
vartable = [x1, x2];

% =============================================
% First, initialize the sum of squares program
prog = sosprogram(vartable);   % No decision variables.

% =============================================
% Next, define the inequality

% p(x1,x2) >=  0
p = 2*x1^4 + 2*x1^3*x2 - x1^2*x2^2 + 5*x2^4;
prog = sosineq(prog,p);

% =============================================
% And call solver
solver_opt.solver = 'sedumi';
[prog,info] = sossolve(prog,solver_opt);

% =============================================
% If program is feasible, p(x1,x2) is an SOS.
echo off;