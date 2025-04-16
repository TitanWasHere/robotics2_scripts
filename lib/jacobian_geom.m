function [J_geom] = jacobian_geom(DH_table, joints, types)
    % INPUT:
    %   - DH_table: table dh, matrix with each row as [theta, alpha, a, d]
    %   - joints: symbolic joint variables, e.g., [q1, q2, ...]
    %   - types: joint types, either 'r' (revolute) or 'p' (prismatic)
    
    % Ensure that the number of joints matches the table and types

    %syms q1 q2 q3 a1 a2 a3 d1
    %DH_table = [
    %    q1, pi/2, a1, d1;
    %    q2, 0, a2,0;
    %    q3, pi/2, a3, 0
    %];
    %types = ["r", "r", "r"];
    %joints = [q1,q2,q3];


    n_rows = size(DH_table, 1);
    n_joints = length(joints);
    if n_rows ~= n_joints
        error("[ERROR]: Number of joints does not match the number of rows in DH_table");
    end
    

    % Initialize the Jacobian matrix
    J_geom = sym(zeros(6, n_joints));
    
    % Initial position and transformation
    T = eye(4);  % Initialize transformation matrix
    mat = DH_HTM(DH_table, "r");
    z_0 = [0; 0; 1];  % The z-axis of the base frame
    p_0_e = mat(1:3, 4);  % End effector position (initially at origin)
    
    % Iterate over each joint to calculate its contribution to the Jacobian
    for i = 1:n_joints
        % Get the transformation matrix for the current joint
       
        % Get the axis of the current joint (z_i-1)
        if i == 1
            z_i_1 = z_0;
        else
            z_i_1 = T(1:3, 1:3) * z_0;  % The z-axis of the previous transformation

        end
        % Create the transformation matrix for this joint
        T_i = DH_HTM(DH_table(i , :), "r");
        T = T * T_i;  % Update the total transformation matrix
        
        
        
        % Determine if the joint is prismatic or revolute
        if types(i) == "r"  % Revolute joint
            % The angular part is z_i-1, the linear part is the derivative of p_0_e
            J_geom(1:3, i) = diff(p_0_e, joints(i));
            J_geom(4:6, i) = z_i_1;
        elseif types(i) == "p"  % Prismatic joint
            % The angular part is zero (0,0,0), the linear part is z_i-1
            J_geom(1:3, i) = z_i_1;
            J_geom(4:6, i) = [0; 0; 0];
        else
            error('Joint type must be "r" (revolute) or "p" (prismatic)');
        end
    end
    J_geom = simplify(J_geom);
end
