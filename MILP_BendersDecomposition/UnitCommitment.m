% Unit Commitment (UC) Example using Benders Decomposition
% Data mostly from [1]
% [1] https://yalmip.github.io/example/unitcommitment/
% [2] G. Morales-Espa?a, J. M. Latorre and A. Ramos, "Tight and Compact MILP Formulation for the Thermal Unit Commitment Problem," in IEEE Transactions on Power Systems, vol. 28, no. 4, pp. 4897-4908, Nov. 2013, doi: 10.1109/TPWRS.2013.22514
clear
clc

%% Basic Conditions
Nunits = 3; % Number of units
Hours = 48; % Number of hours

Pmax = [100; 50; 25]; % Maximum power output for each unit
Pmin = [20; 40; 1]; % Minimum power output for each unit

C = [10 20 20]; % Cost per unit power for each unit
Cu_NL = [2 3 4]; % No-load cost for each unit
Cup = [5 10 5]; % Startup cost for each unit
Cdown = [5 10 5]; % Shutdown cost for each unit

Pforecast = 100 + 50*sin((1:Hours)*2*pi/24); % Power forecast for each hour

%% Variable Definition
u = binvar(Nunits, Hours, 'full'); % Binary variable indicating if unit g is on at time t
v = binvar(Nunits, Hours, 'full'); % Binary variable indicating if unit g starts up at time t
w = binvar(Nunits, Hours, 'full'); % Binary variable indicating if unit g shuts down at time t
P = sdpvar(Nunits, Hours, 'full'); % Power output of unit g at time t

%% Constraints
Cons = [];

% Constraint (1): Power output limits when unit is on
for t = 1:Hours
    Cons = [Cons, u(:,t).*Pmin <= P(:,t) <= u(:,t).*Pmax];
end

% Adding minimum up- and down-time constraints
TUg = [6; 30; 1]; % Minimum uptime for each unit
TDg = [3; 6; 3]; % Minimum downtime for each unit

% Constraint (6) in [2]: Ensure minimum uptime
for g = 1:Nunits
    for t = TUg(g):Hours
        Cons = [Cons, sum(v(g, (t-TUg(g)+1):t), 2) <= u(g, t)];
    end
end

% Constraint (7) in [2]: Ensure minimum downtime
for g = 1:Nunits
    for t = TDg(g):Hours
        Cons = [Cons, sum(w(g, (t-TDg(g)+1):t), 2) <= 1 - u(g, t)];
    end
end

% Constraint (8) in [2]: Link startup and shutdown variables with unit status
for t = 2:Hours
    Cons = [Cons, u(:,t) - u(:,t-1) == v(:,t) - w(:,t)];
end

% Constraint (9): Meet power demand
for t = 1:Hours
    Cons = [Cons, sum(P(:,t)) >= Pforecast(t)];
end

% Redundant constraint to limit the "complicating variables" y
for t = 1:Hours
    Cons = [Cons, sum(u(:,t).*Pmax) >= Pforecast(t)];
end

%% Objective Function
Obj_f = 0; % Fuel cost
Obj_u = 0; % No-load cost
Obj_up = 0; % Startup cost
Obj_down = 0; % Shutdown cost
for t = 1:Hours
    Obj_f = Obj_f + C*P(:,t);
    Obj_u = Obj_u + Cu_NL*u(:,t);
    Obj_up = Obj_up + Cup*v(:,t);
    Obj_down = Obj_down + Cdown*w(:,t);
end
Obj = Obj_f + Obj_u + Obj_up + Obj_down;

%% (1) Solve the UC problem by Gurobi
ops = sdpsettings('solver','gurobi','verbose',0);
result_NoBD = optimize(Cons, Obj, ops);
s_P = value(P);
s_u = value(u);
s_v = value(v);
s_w = value(w);
s_Obj = value(Obj);

figure
h1 = bar(s_P', 'stack');
legend('Unit 1', 'Unit 2', 'Unit 3');
title('UC Result Solved by Gurobi');

disp(['Solved by Gurobi, without Benders Decomposition: ', num2str(result_NoBD.solvertime), ' s']);

%% (2) Solve the Problem by Benders Decomposition
% 2.1 Export Model with Gurobi
ops = sdpsettings('solver','gurobi','verbose',1);
[model, r_model, diagnostic, internalmodel] = export(Cons, Obj, ops);% 导出模型的矩阵参数(模型是Yamlip形式)

% 2.2 Solve it by BendersDecomposition (embedded with Gurobi)
BendersDecomposition;