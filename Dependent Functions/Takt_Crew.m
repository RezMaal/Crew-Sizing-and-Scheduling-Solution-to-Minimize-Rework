% Constants
function [result,model]=Takt_Crew(k,n,C_c,C_o,C_t,a,M_col,M_r)

% Initialize Gurobi model
model = struct();
model.modelsense = 'min'; % Set as a minimization problem
model.obj = []; % Initialize the objective coefficients
model.A = sparse([]); % Matrix to hold constraints
model.rhs = []; % Right-hand side of constraints
model.sense = ''; % Inequality/equality sense for constraints

% Variable definitions
num_X_vars = k * n; % Total number of X variables
num_b_vars = k * n; % Total number of b variables
num_Z_vars = k; % Total number of Z variables
num_Y_vars = k * k * n; % Total number of Y variables

total_vars = num_X_vars + num_b_vars + num_Z_vars + num_Y_vars;

% Initialize variables
Mx=max(max(M_col),max(M_r));
model.lb = zeros(total_vars, 1); % Lower bounds (0 for all variables)
model.vtype = repmat('I', total_vars, 1); % Variable types ('C' for continuous, 'I' for integer, 'B' for binary)

% Define the variable indices
X_idx = 1:num_X_vars;
b_idx = (1:num_b_vars) + num_X_vars;
Z_idx = (1:num_Z_vars) + num_X_vars + num_b_vars;
Y_idx = (1:num_Y_vars) + num_X_vars + num_b_vars + num_Z_vars;

% Adjust variable types
model.vtype(X_idx) = 'I'; % X is integer
model.vtype(b_idx) = 'B'; % b is binary
model.vtype(Z_idx) = 'B'; % Z is binary
model.vtype(Y_idx) = 'I'; % X is integer
model.ub(X_idx) = Mx*ones(num_X_vars, 1); % Upper bounds (Mx for integer, and 1 for binary)
model.ub(b_idx) = ones(num_b_vars, 1);
model.ub(Y_idx) = Mx*ones(num_Y_vars, 1); % Upper bounds (Mx for integer, and 1 for binary)
model.ub(Z_idx) = ones(num_Z_vars, 1);

% Objective function
obj_X = C_c * ones(num_X_vars, 1);
obj_b = -C_t * ones(num_b_vars, 1)/n;
obj_Z = C_o * ones(num_Z_vars, 1) + C_t * ones(num_Z_vars, 1);
obj_Y = zeros(num_Y_vars, 1); % Y does not appear in the objective directly

model.obj = [obj_X; obj_b; obj_Z; obj_Y];

% Constraints
row = 0;

% Add constraints for X, b, and Z
for i = 1:k
    % Row-specific M_r(i)
    M_r_i = M_r(i);
    
    % Row indices for X and b
    ix=sub2ind([k n],repmat(i,1,n),(1:n));
    X_row_idx = X_idx(ix);
    b_row_idx = b_idx(ix);
    
    % Constraint: Z(i) based on sum of row of X
    row = row + 1;
    model.A(row, X_row_idx) = 1; % Sum of X(i, :)
    model.A(row, Z_idx(i)) = -M_r_i; % -M_r(i) * Z(i)
    model.rhs(row) = 0;
    model.sense(row) = '<';
    
    row = row + 1;
    model.A(row, X_row_idx) = 1; % Sum of X(i, :)
    model.A(row, Z_idx(i)) = -1; % Z(i)
    model.rhs(row) = 0;
    model.sense(row) = '>';
    
    % Constraint: Z(i) based on sum of row of b
    row = row + 1;
    model.A(row, b_row_idx) = 1; % Sum of b(i, :)
    model.A(row, Z_idx(i)) = -1; % Z(i)
    model.rhs(row) = 0;
    model.sense(row) = '>';
    
    % Add constraints for b with respect to X
    for j = 1:n
        % Indices for X(i, j) and b(i, j)
        ixb=sub2ind([k n],i,j);
        X_ij_idx = X_idx(ixb);
        b_ij_idx = b_idx(ixb);
        
        % Constrain b based on X
        row = row + 1;
        model.A(row, b_ij_idx) = 1; % b(i, j)
        model.A(row, X_ij_idx) = -1; % -X(i, j)
        model.rhs(row) = 0;
        model.sense(row) = '<'; % b(i, j) <= X(i, j)
        
        row = row + 1;
        model.A(row, b_ij_idx) = -M_col(j); % -b(i, j)
        model.A(row, X_ij_idx) = 1 ; % X(i, j) / M_col(j)
        model.rhs(row) = 0;
        model.sense(row) = '<'; % b(i, j) >= X(i, j) / M_col(j)
    end
end

% Add constraints for Y
for j = 1:n
    M_j = M_col(j); % Column-specific big-M
    
    % Sum constraint for Y_j
    iyj=sub2ind([k k n],1:k,1:k,repmat(j,1,k));
    Y_col_idx = Y_idx(iyj); % Indices for all elements in Y_j
    row = row + 1;
    model.A(row, Y_col_idx) = 1; % Sum of all elements in Y_j
    model.rhs(row) = a(j); % Upper bound is a_j
    model.sense(row) = '>';
    
    for i = 1:k
        for l = 1:k
            % Indices for Y(i, l, j)
            ix=sub2ind([k n],i,j);
            ib=sub2ind([k n],l,j);
            iy=sub2ind([k k n],i,l,j);
            % Y_ij_idx = Y_idx((j - 1) * k * k + (i - 1) * k + l);
            Y_ij_idx = Y_idx(iy);
            X_ij_idx = X_idx(ix);
            b_lj_idx = b_idx(ib);
            
            % Add linearization constraints for Y(i, l, j)
            row = row + 1;
            model.A(row, Y_ij_idx) = 1; % Y(i, l, j)
            model.A(row, X_ij_idx) = -1; % -X(i, j)
            model.rhs(row) = 0;
            model.sense(row) = '<';
            
            row = row + 1;
            model.A(row, Y_ij_idx) = 1; % Y(i, l, j)
            model.A(row, b_lj_idx) = -M_j; % -M_j * b(l, j)
            model.rhs(row) = 0;
            model.sense(row) = '<';
            
            row = row + 1;
            model.A(row, Y_ij_idx) = -1; % -Y(i, l, j)
            model.A(row, X_ij_idx) = 1; % X(i, j)
            model.A(row, b_lj_idx) = M_j; % M_j * (1 - b(l, j))
            model.rhs(row) = M_j;
            model.sense(row) = '<';

            row = row + 1;
            model.A(row, Y_ij_idx) = -1; % -Y(i, l, j)
            model.A(row, b_ij_idx) = 1; % b(l, j)
            model.rhs(row) = 0;
            model.sense(row) = '<';
        end
    end
end

% Solve the model using Gurobi
params.outputflag = 1; % Enable solver output
params.FeasibilityTol = 1e-9;
params.IntFeasTol     = 1e-9;
params.OptimalityTol  = 1e-9;
params.MIPGap         = 1e-8;
result = gurobi(model, params);

% Display results
if strcmp(result.status, 'OPTIMAL')
    fprintf('Optimal objective value: %f\n', result.objval);
else
    fprintf('No optimal solution found. Status: %s\n', result.status);
end