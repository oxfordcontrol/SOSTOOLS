% SOSDEMO7s --- Chebyshev polynomials
% Section 4.7 of SOSTOOLS User's Manual

clear; echo on;
syms x gam;

% Degree of Chebyshev polynomial
ndeg = 8;   

% =============================================
% First, initialize the sum of squares program
prog = sosprogram([x],[gam]);

% Create the polynomial P
Z = monomials(x,[0:ndeg-1]);
[prog,P1] = sospolyvar(prog,Z);
P = P1 + gam * x^ndeg;           % The leading coeff of P is gam

% Imposing the inequalities
prog = sosineq(prog, 1 - P, [-1, 1]);
prog = sosineq(prog, 1 + P, [-1, 1]);

% And setting objective
prog = sossetobj(prog, -gam);

% Then solve the program
solver_opt.solver = 'sedumi';
prog = sossolve(prog,solver_opt);

% =============================================
% Finally, get solution
SOLV = sosgetsol(prog,P)
GAM = sosgetsol(prog,gam)
echo off