% Description: Implement McCormick envelope linearization using Yalmip and
% Gurobi solver, a Simple Example.

% Author: X Yang 
% Date: 2025/4/17
% E-mail: 
% Ref. : https://blog.csdn.net/qq_40926887/article/details/114789742

clear all;
close all;
clc;

% Check if YALMIP is installed
if exist('sdpvar') ~= 2
    error('Please install YALMIP toolbox first');
end

% Define decision variables
x = sdpvar(1);
y = sdpvar(1);
w = sdpvar(1); % Represents the linearized term x*y

% Set variable bounds
x_lb = -2; x_ub = 2;   % x ¡Ê [-2, 2]
y_lb = -1; y_ub = 3;   % y ¡Ê [-1, 3]

% McCormick envelope constraints in matrix form
% Reformulate inequalities as A*[w;x;y] <= b
A = [-1,  y_lb,  x_lb;   % -w + y_lb*x + x_lb*y <= x_lb*y_lb (Lower bound 1)
     -1,  y_ub,  x_ub;   % -w + y_ub*x + x_ub*y <= x_ub*y_ub (Lower bound 2)
      1, -y_lb, -x_ub;   % w - y_lb*x - x_ub*y <= -x_ub*y_lb (Upper bound 1)
      1, -y_ub, -x_lb];  % w - y_ub*x - x_lb*y <= -x_lb*y_ub (Upper bound 2)

b = [x_lb*y_lb;
     x_ub*y_ub;
     -x_ub*y_lb;
     -x_lb*y_ub];

% Construct constraints
constraints = [A*[w;x;y] <= b];
constraints = [constraints, 
               x_lb <= x <= x_ub, 
               y_lb <= y <= y_ub];

% Define objective function (minimize w)
objective = w;

% Configure Gurobi solver
ops = sdpsettings('solver', 'gurobi', 'verbose', 1);

% Solve optimization problem
result = optimize(constraints, objective, ops);

% Display results
if result.problem == 0
    disp('Optimal solution found:');
    disp(['w = ', num2str(value(w))]);
    disp(['x = ', num2str(value(x))]);
    disp(['y = ', num2str(value(y))]);
else
    disp('Optimization failed');
end