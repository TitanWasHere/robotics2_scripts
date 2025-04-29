clc; clear; close all;
syms alpha a d theta real;

%% Define symbolic DH parameters (General Form)
paramOrder = {'alpha', 'a', 'd', 'theta'};

% Define the D-H Table Symbolically (Modify this to define any robot)
DH_sym = [
    sym('alpha1'), sym('a1'), sym('d1'), sym('theta1');
    sym('alpha2'), sym('a2'), sym('d2'), sym('theta2');
    sym('alpha3'), sym('a3'), sym('d3'), sym('theta3')
];

% Specify the uncertain parameters that need to be estimated
uncertain_params = ["a1", "d1", "theta2", "a3"]; % Modify as needed

% Compute Direct Kinematics Symbolically
T_sym = directKinematics(DH_sym, paramOrder);
disp('Symbolic Transformation Matrix:');
disp(T_sym);

% Extract symbolic end-effector position
ee_sym = T_sym(1:3, 4);
disp('Symbolic End-Effector Position:');
disp(ee_sym);

%% Convert to Numeric Values for Calibration
n_joints = size(DH_sym, 1);

% Define nominal values (Initial guess)
DH_nom = [ % Change these values according to the robot's nominal parameters
    pi/2, 0.3, 0.2, pi/4;
    0, 0.5, 0, pi/3;
    -pi/2, 0.2, 0, pi/6
];

% Define true values (introducing small manufacturing errors)
DH_true = DH_nom + 0.01 * randn(size(DH_nom));

% Generate calibration samples dynamically
num_samples = 20;
joint_types = [0, 1, 1];
joint_samples = generateJointSamples(DH_true, joint_types, num_samples);

% Compute ground truth using true DH parameters
ee_true = zeros(3, num_samples);
for i = 1:num_samples
    ee_true(:, i) = double(subs(ee_sym, DH_sym(:), DH_true(:)));
end

%% Calibration Algorithm (Iterative Least Squares)
e_r = 1e-4;
i_max = 50;
mse_history = zeros(i_max, 1);
DH_cal = DH_nom;

% Compute symbolic regressor matrix (Jacobian of EE position w.r.t. only uncertain parameters)
uncertain_syms = sym(uncertain_params);
Phi_sym = jacobian(ee_sym, uncertain_syms);

for iter = 1:i_max
    % Compute estimated end-effector positions with current DH parameters
    ee_est = zeros(3, num_samples);
    for i = 1:num_samples
        ee_est(:, i) = double(subs(ee_sym, DH_sym(:), DH_cal(:)));
    end

    % Compute error
    delta_r = ee_true - ee_est;
    
    % Compute Mean Squared Error (MSE)
    mse = mean(delta_r(:).^2);
    mse_history(iter) = mse;
    fprintf('Iteration %d, MSE: %.6f\n', iter, mse);
    
    % Stopping criterion
    if mse < e_r
        fprintf('Converged in %d iterations.\n', iter);
        break;
    end

    % Evaluate numerical regressor matrix with only uncertain parameters
    Phi = zeros(3*num_samples, numel(uncertain_syms));
    for i = 1:num_samples
        Phi(3*i-2:3*i, :) = double(subs(Phi_sym, DH_sym(:), DH_cal(:)));
    end

    % Solve for parameter corrections using least squares
    delta_phi = pinv(Phi) * delta_r(:);
    
    % Update only the uncertain parameters in DH_cal
    for j = 1:numel(uncertain_syms)
        [row, col] = find(DH_sym == uncertain_syms(j));
        DH_cal(row, col) = DH_cal(row, col) + delta_phi(j);
    end
end

% Display results
disp('Calibrated D-H Parameters:');
disp(DH_cal);

% Plot MSE convergence
figure;
plot(1:iter, mse_history(1:iter), '-o');
xlabel('Iteration');
ylabel('Mean Squared Error (MSE)');
title('Convergence of Calibration');
grid on;

%% Function: Direct Kinematics
function T = directKinematics(DH, paramOrder)
    n = size(DH, 1);
    paramMap = containers.Map({'alpha', 'a', 'd', 'theta'}, 1:4);
    paramIdx = cellfun(@(x) paramMap(x), paramOrder);
    T = eye(4);

    for i = 1:n
        alpha = DH(i, paramIdx(1));
        a = DH(i, paramIdx(2));
        d = DH(i, paramIdx(3));
        theta = DH(i, paramIdx(4));

        Ti = [cos(theta), -sin(theta)*cos(alpha), sin(theta)*sin(alpha), a*cos(theta);
              sin(theta), cos(theta)*cos(alpha), -cos(theta)*sin(alpha), a*sin(theta);
              0, sin(alpha), cos(alpha), d;
              0, 0, 0, 1];

        T = T * Ti;
        T = simplify(T);
    end
end

%% Function: Generate Joint Samples
function jointSamples = generateJointSamples(DH, jointTypes, numSamples)
    jointSamples = zeros(size(DH, 1), numSamples);
    
    for i = 1:size(DH, 1)
        if jointTypes(i) == 1  % Revolute
            jointSamples(i, :) = (rand(1, numSamples) - 0.5) * pi;
        else  % Prismatic
            nominal_d = DH(i, 3);
            jointSamples(i, :) = nominal_d + (rand(1, numSamples) - 0.5) * 0.1;
        end
    end
end