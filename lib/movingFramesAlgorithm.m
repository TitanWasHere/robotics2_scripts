function [v, w, vc] = movingFramesAlgorithm(DH, dq, rc, joint_types)
% movingFramesAlgorithm computes the velocity and acceleration of a
% manipulator using the moving frames algorithm.
% Inputs:
%   DH: Denavit-Hartenberg parameters (Nx4 matrix) [alpha, a, d, theta]
%   qd: symbolic joint velocities (Nx1 vector)
%   rc: position of the center of mass of each link(Nx3 vector)
%   joint_types: types of joints (Nx1 vector) ["r": revolute, "p": prismatic]
% Outputs:
%   v: linear velocity of the end-effector (3x1 vector)
%   w: angular velocity of the end-effector (3x1 vector)
%   vc: velocity of the center of mass of each link (3xN matrix)

    n = size(DH, 1);
    v = cell(n, 1); 
    w = cell(n, 1); 
    T = cell(n, 1);
    A = cell(n, 1);
    R = cell(n, 1); 
    vc = cell(n,1);
    r = cell(n,1);
    p = cell(n,1);
    sigma = cell(n,1);
    DK = cell(n,1);

    w_prev = zeros(3, 1);
    v_prev = zeros(3, 1);
    
    T_prev = eye(4);

    for i = 1:n
        alpha = DH(i, 1);
        a = DH(i, 2);
        d = DH(i, 3);
        theta = DH(i, 4);

        % Compute the transformation matrix
        T_i = [cos(theta),    -sin(theta)*cos(alpha),     sin(theta)*sin(alpha),  a*cos(theta);
             sin(theta),    cos(theta)*cos(alpha),      -cos(theta)*sin(alpha), a*sin(theta);
             0,             sin(alpha),                 cos(alpha),             d;
             0,             0,                          0,                      1];
        A{i} = T_i;
        DK{i} = T_prev * T_i;
        R{i} = A{i}(1:3, 1:3);
        p{i} = A{i}(1:3, 4);
        
        r{i} = R{i}' * p{i};
        
        T_prev = T_i;

        if joint_types(i) == "r" || joint_types(i) == "R"
            sigma{i} = 0;
        elseif joint_types(i) == "p" || joint_types(i) == "P"
            sigma{i} = 1;
        else
            error("Invalid joint type. Use 'r' for revolute and 'p' for prismatic.");
        end

        z0 = [0; 0; 1];
        w{i} = simplify(R{i}' * (w_prev + (1 - sigma{i}) * dq(i) * z0));
        w_prev = w{i};

        v{i} = simplify(R{i}' * (v_prev + sigma{i} * dq(i) * z0) + cross(w{i}, r{i}));
        v_prev = v{i};

        vc{i} = simplify(v{i} + cross(w{i}, rc(:, i)));
        
    end

    for i = 1:n
        fprintf("Link %d:\n", i);
        fprintf("Angular velocity (w):\n"); disp(w{i});
        fprintf("Linear velocity (v):\n"); disp(v{i});
        fprintf("Velocity of center of mass (vc):\n"); disp(vc{i});
        fprintf("Rotation matrix (R):\n"); disp(R{i});
        fprintf("point (r):\n"); disp(r{i});
        fprintf("cross product (w x r):\n"); disp(cross(w{i}, r{i}));
        fprintf("Sigma: %d\n", sigma{i});
        fprintf("\n")
    end

end