function cubicTrajectory(q_i, q_f, t_i, t_f, v_i, v_f)
    % cubicTrajectory - Computes and plots a cubic trajectory for multiple robot joints.
    %
    % Syntax: cubicTrajectory(q_i, q_f, t_i, t_f, v_i, v_f)
    %
    % Inputs:
    %   q_i  - Initial joint positions (vector)
    %   q_f  - Final joint positions (vector)
    %   t_i  - Initial time (seconds)
    %   t_f  - Final time (seconds)
    %   v_i  - Initial joint velocities (vector)
    %   v_f  - Final joint velocities (vector)
    %
    % Outputs:
    %   Plots of position, velocity, and acceleration profiles for each joint.
    %   Displays the maximum velocity for each joint.
    %
    % Example:
    %   cubicTrajectory([0, 0], [1, 1], 0, 2, [0, 0], [0.5, 0.5]);

    %% Parameters
    T = t_f - t_i; % Total duration
    Delta_q = q_f - q_i;
    numJoints = length(q_i);

    %% Determine the cubic coefficients
    a0 = q_i;
    a1 = v_i;
    a3 = (v_f - v_i - 2 .* Delta_q / T) ./ T^2;
    a2 = (3 .* Delta_q / T - 2 .* v_i - v_f) ./ T;

    %% Define the symbolic variable and expressions
    syms t real
    tau = t - t_i;  % time shift so that tau=0 corresponds to t_i

    % Position, Velocity, and Acceleration expressions for each joint
    q_sym = a0 + a1*tau + a2*tau^2 + a3*tau^3;
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
    q_vals   = cell2mat(q_vals');
    qd_vals  = cell2mat(qd_vals');
    qdd_vals = cell2mat(qdd_vals');


    %% Compute and display maximum velocity and corresponding time for each joint
    max_velocity = zeros(1, numJoints);
    time_of_max_velocity = zeros(1, numJoints);
    for j = 1:numJoints
        [max_velocity(j), idx] = max(abs(qd_vals(:, j)));
        time_of_max_velocity(j) = time_vec(idx);
    end
    disp('Maximum velocity for each joint:');
    disp(max_velocity);
    disp('Time at which maximum velocity is reached for each joint:');
    disp(time_of_max_velocity);

    %% Plot the results
    figure;
    for j = 1:numJoints
        subplot(3, numJoints, j);
        plot(time_vec, q_vals(:, j), 'b', 'LineWidth', 2);
        grid on;
        xlabel('Time (s)');
        ylabel(['Position Joint ', num2str(j)]);
        title(['Joint ', num2str(j), ' Position Profile']);
        set(gca, 'FontSize', 12);

        subplot(3, numJoints, j + numJoints);
        plot(time_vec, qd_vals(:, j), 'r', 'LineWidth', 2);
        grid on;
        xlabel('Time (s)');
        ylabel(['Velocity Joint ', num2str(j)]);
        title(['Joint ', num2str(j), ' Velocity Profile']);
        set(gca, 'FontSize', 12);

        subplot(3, numJoints, j + 2*numJoints);
        plot(time_vec, qdd_vals(:, j), 'k', 'LineWidth', 2);
        grid on;
        xlabel('Time (s)');
        ylabel(['Acceleration Joint ', num2str(j)]);
        title(['Joint ', num2str(j), ' Acceleration Profile']);
        set(gca, 'FontSize', 12);
    end

    set(gcf, 'Position', [100, 100, 1200, 800]); % Increase figure width for multiple joints

    %% Display the evaluated expressions (optional)
    disp('Position q(t):');
    disp(vpa(q_sym, 6));
    disp('Velocity qdot(t):');
    disp(vpa(qd_sym, 6));
    disp('Acceleration qddot(t):');
    disp(vpa(qdd_sym, 6));
end
