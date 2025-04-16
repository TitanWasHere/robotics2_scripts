function quinticTrajectory(q_i, q_f, t_i, t_f, v_i, v_f, a_i, a_f)
    % quinticTrajectory - Computes and plots a quintic trajectory for multiple robot joints.
    %
    % Syntax: quinticTrajectory(q_i, q_f, t_i, t_f, v_i, v_f, a_i, a_f)
    %
    % Inputs:
    %   q_i  - Initial joint positions (vector)
    %   q_f  - Final joint positions (vector)
    %   t_i  - Initial time (seconds)
    %   t_f  - Final time (seconds)
    %   v_i  - Initial joint velocities (vector)
    %   v_f  - Final joint velocities (vector)
    %   a_i  - Initial joint accelerations (vector)
    %   a_f  - Final joint accelerations (vector)
    %
    % Example:
    %   quinticTrajectory([-pi/4 ; pi/4 ; pi/4], [0;0;pi/4], 0, 2, [1;-1;0], [0;0;0], [0;0;0], [0;0;0]);

    %% Parameters
    T = t_f - t_i; % Total duration
    Delta_q = q_f - q_i;
    numJoints = length(q_i);

    %% Define the symbolic variable and expressions
    syms t real
    tau = t - t_i;  % time shift so that tau=0 corresponds to t_i

    % Here, we use a formulation that blends the initial and final conditions
    % into a smooth trajectory.
    q_sym = (1-tau/T)^3 * (q_i + (3 * q_i + v_i * T) * (tau/T) + ((a_i * T^2 + 6 * v_i * T + 12 * q_i)/2) * (tau/T)^2) ...
          + (tau/T)^3 * (q_f + (3*q_f - v_f * T)*(1-tau/T) + ((a_f * T^2 - 6 * v_f * T + 12*q_f)/2) * (1-tau/T)^2);
    qd_sym = diff(q_sym, t);
    qdd_sym = diff(qd_sym, t);

    %% Convert symbolic expressions to MATLAB functions for evaluation
    q_fun   = matlabFunction(q_sym, 'Vars', t);
    qd_fun  = matlabFunction(qd_sym, 'Vars', t);
    qdd_fun = matlabFunction(qdd_sym, 'Vars', t);

    %% Time vector for plotting
    numPoints = 100;
    time_vec = linspace(t_i, t_f, numPoints);

    %% Evaluate the profiles
    q_vals   = arrayfun(q_fun, time_vec, 'UniformOutput', false);
    qd_vals  = arrayfun(qd_fun, time_vec, 'UniformOutput', false);
    qdd_vals = arrayfun(qdd_fun, time_vec, 'UniformOutput', false);
    
    % Convert the cell arrays to matrices and transpose so that each column is a joint's data
    q_vals   = cell2mat(q_vals)';    % Now size is 100xnumJoints
    qd_vals  = cell2mat(qd_vals)';     % Now size is 100xnumJoints
    qdd_vals = cell2mat(qdd_vals)';    % Now size is 100xnumJoints

    %% Compute and display maximum velocity for each joint
    max_velocity = max(abs(qd_vals));
    disp('Maximum velocity for each joint:');
    disp(max_velocity);

    %% Plot the results in a single row (horizontally arranged) and square each subplot
    figure;
    joint_colors = ['b', 'r', 'k']; % Define colors for joints (adjust as needed)

    % Position plot
    subplot(1,3,1);
    hold on;
    for j = 1:numJoints
        plot(time_vec, q_vals(:, j), joint_colors(j), 'LineWidth', 2);
    end
    grid on;
    xlabel('Time [s]');
    ylabel('[rad]');
    title('Joint Positions');
    legend('Joint 1', 'Joint 2', 'Joint 3');
    axis square;   % Forces the subplot to have equal axis lengths
    set(gca, 'FontSize', 12);

    % Velocity plot
    subplot(1,3,2);
    hold on;
    for j = 1:numJoints
        plot(time_vec, qd_vals(:, j), joint_colors(j), 'LineWidth', 2);
    end
    grid on;
    xlabel('Time [s]');
    ylabel('[rad/s]');
    title('Joint Velocities');
    legend('Joint 1', 'Joint 2', 'Joint 3');
    axis square;
    set(gca, 'FontSize', 12);

    % Acceleration plot
    subplot(1,3,3);
    hold on;
    for j = 1:numJoints
        plot(time_vec, qdd_vals(:, j), joint_colors(j), 'LineWidth', 2);
    end
    grid on;
    xlabel('Time [s]');
    ylabel('[rad/s²]');
    title('Joint Accelerations');
    legend('Joint 1', 'Joint 2', 'Joint 3');
    axis square;
    set(gca, 'FontSize', 12);

    % Adjust figure size (square-ish overall layout)
    set(gcf, 'Position', [100, 100, 1200, 400]); 

    %% Display the evaluated expressions (optional)
    disp('Position q(t):');
    disp(vpa(simplify(q_sym), 4));
    disp('Velocity qdot(t):');
    disp(vpa(simplify(qd_sym), 4));
    disp('Acceleration qddot(t):');
    disp(vpa(simplify(qdd_sym), 4));
end
