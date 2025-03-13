%% Benders Decomposition Algorithm for Mixed Integer Linear Programming (MILP)
% Author: Xiang Yang
% Date: 2025-03-13
% Reference: 
% [1] https://github.com/a280558071/UC_Solvedby_BendersDecomp/blob/main/BendersDecomposition_Gurobi_2022.m

%% Standard Form of the Original Problem:
% min d'*y + c'*x
% s.t. F*y + E*x <= h;
%      A*y <= b
%      lx <= x <= ux
%      ly <= y <= uy
%      y: binaries (complicating variables)
%      x(>=0): continuous

%% Master Problem (MP1):
% min z
% s.t. z >= d'*y
%      A*y <= b
%      y: binaryies

%% Subproblem (SP1) for (1-y_hat):
% min c'*x
% s.t. E*x <= h - F*y_hat
%      x >= 0

%% Dual of Subproblem (SP2):
% max (F*y_hat)'*u
% s.t. E'*u <= c
%      u <= 0
t_BD_s = tic;  
% Import Gurobi model and define parameters
epsilon = 1e-5; % Convergence tolerance
z_UB = inf; % Upper bound of the objective value
z_LB = -inf; % Lower bound of the objective value
MipDisplayInterval = 10000; % Display interval for MIP solving output
MaxIter = 100; % Maximum number of iterations
Outputflag = 0; % 0 means no MIP solving output

% Transform variables into Benders decomposition form
Ind_y = []; % Indices of binary variables (complicating variables)
Ind_x = []; % Indices of continuous variables
for i = 1:length(model.vtype)
    if strcmp(model.vtype(i), 'B')
        Ind_y = [Ind_y, i];
    end
    if strcmp(model.vtype(i), 'C')
        Ind_x = [Ind_x, i];
    end
end
N_y = length(Ind_y);
N_x = length(Ind_x);
N_p = N_y + N_x; % Total number of variables in the primal problem
N_d = length(model.sense); % Number of constraints in the dual problem

% Extract coefficients and right-hand sides from the model
FA = model.A(:, 1:N_y); % Coefficients of binary variables in A
E = model.A(:, (N_y+1):end); % Coefficients of continuous variables in A
hb = model.rhs;
d = model.obj(Ind_y);
c = model.obj(Ind_x);

% Separate constraints into master problem (MP) and subproblem (SP)
MP_rows = [];
SP_rows = [];
for i = 1:N_d
    if isempty(find(E(i,:)))
        MP_rows = [MP_rows, i];
    else
        SP_rows = [SP_rows, i];
    end
end
A = FA(MP_rows, :);
F = FA(SP_rows, :);
b = hb(MP_rows);
h = hb(SP_rows);

% Step 1: Solve Master Problem (MP1)
MP.obj = [1; zeros(N_y, 1)]; % Objective function for MP1: minimize z subject to constraints
MP.A = sparse([1, -d']); % Constraints: z >= d'y (Optimality cut constraint)
MP.A = [MP.A; zeros(length(MP_rows), 1), A]; % Add remaining constraints from A*y <= b (Original constraints involving y only)
MP.sense = ['>'; model.sense(MP_rows)]; % Constraint senses: '>' for optimality cut, others as in original model
MP.rhs = [0; b]; % Right-hand sides
MP.lb = [0; model.lb(Ind_y)]; % Variable lower bounds: z >= 0, y within their original bounds
MP.ub = [Inf; model.ub(Ind_y)]; % Variable upper bounds: z unbounded, y within their original bounds
MP.vtype = ['C'; model.vtype(Ind_y)]; % Variable types: z continuous, y binary
MP.params.DisplayInterval = MipDisplayInterval; % Set display interval
MP.params.outputflag = Outputflag; % Suppress solver output
r_MP = gurobi(MP, MP.params); % Solve MP1
z_LB = r_MP.objval; % Update lower bound

% Initialize counters for cuts
p = 0;
q = 0;

% Construct Subproblem (SP1) or its dual (SP2)
SP.A = E(SP_rows, :);
SP.obj = c;
SP.sense = model.sense(SP_rows);
SP.lb = model.lb(Ind_x);
SP.ub = model.ub(Ind_x);
SP.vtype = model.vtype(Ind_x);
SP.params.DisplayInterval = MipDisplayInterval;
SP.params.outputflag = Outputflag;
r_SP = gurobi(SP, SP.params); % Solve SP1 initially

timeS = 0;
iter = 1;
s_P_BD = [];
s_u_BD = [];
s_v_BD = [];
s_w_BD = [];
abs_error = abs((z_UB - z_LB) / z_UB);

% Main loop of Benders decomposition
while iter <= MaxIter
    % Step 2: Solve Subproblem (SP1) with current y_hat
    y_hat = r_MP.x(2:end);
    SP.rhs = h - F * y_hat;
    r_SP = gurobi(SP, SP.params);
    
    if strcmp(r_SP.status, 'OPTIMAL')
        % Add optimality cut
        p = p + 1;
        fprintf('Add optimality cut %d!\n', p);
        z_UB_new = r_SP.objval + d' * y_hat;
        if z_UB_new <= z_UB
            z_UB = z_UB_new;
            s_u_BD_final = r_SP.pi;
            s_v_BD_final = r_SP.x;
            s_w_BD_final = r_SP.slack;
            s_P_BD_final = r_SP.x;
        end
        abs_error = abs((z_UB - z_LB) / z_UB);
        fprintf('Upper Bound: %.4f  Lower Bound: %.4f  Gap: %.2f%%\n', z_UB, z_LB, round(abs_error * 100, 2));
        if abs_error <= epsilon
            fprintf('Benders decomposition converged! The final results are:\n');
            fprintf('Upper Bound: %.4f  Lower Bound: %.4f  Gap: %.2f%%\n', z_UB, z_LB, round(abs_error * 100, 2));
            break
        end
        DualVar = r_SP.pi;
        oc(p).A = [1, (DualVar' * F - d')]; % Optimality cut: z >= d'y + (h - F*y)'*pi -> z + (pi'*F - d')*y >= pi'*h
        oc(p).b = DualVar' * h;
        MP.A = [MP.A; oc(p).A];
        MP.rhs = [MP.rhs; oc(p).b];
        MP.sense = [MP.sense; '>'];
    elseif strcmp(r_SP.status, 'INF_OR_UNBD')
        % Add feasibility cut
        SP_n = SP;
        SP_n.obj = [zeros(N_x, 1); ones(length(SP_rows), 1)];
        SP_n.A = [SP.A, -diag(ones(length(SP_rows), 1))];
        SP_n.lb = [SP.lb; zeros(length(SP_rows), 1)];
        SP_n.ub = [SP.ub; ones(length(SP_rows), 1) * Inf];
        SP_n.vtype = [SP.vtype; ones(length(SP_rows), 1) * 'C'];
        r_SP_n = gurobi(SP_n, SP.params);
        q = q + 1;
        fprintf('Add feasibility cut %d!\n', q);
        DualVar = r_SP_n.pi; % Dual variables u^r
        fc(q).A = [0, DualVar' * F]; % Feasibility cut: (h - F*y)'*u^r <= 0 -> u^r'*F*y >= u^r'*h
        fc(q).b = DualVar' * h;
        MP.A = [MP.A; fc(q).A];
        MP.rhs = [MP.rhs; fc(q).b];
        MP.sense = [MP.sense; '>'];
    end
    
    % Step 3: Solve updated Master Problem (MP2)
    r_MP = gurobi(MP, MP.params);
    z_LB = r_MP.objval;
    iter = iter + 1;
end

t_BD_e = toc(t_BD_s); % End timer
fprintf('Computation Time (s): %.2f s\n', t_BD_e);



