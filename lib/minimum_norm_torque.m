function tau = minimum_norm_torque(J,  J_dot, q_dot, M, r_ddot, n, tau_0)
    %Inputs:
    %   J : (m x n) jacobian matrix of the task 
    %   J_dot : (m x n) time derivative of the jacobian matrix
    %   q_dot : (m x 1) velocity of the joints
    %   M : (n x n) inertia matrix
    %   n :  (n x 1) Coriolis and centrifugal forces, gravity vector, and other external forces
    %   tau_0 : (n x 1) nullspace 
    %Output:
    %   tau : (n x 1) minimum torque norm solution for dynamic redundancy
    
    x_ddot = r_ddot - J_dot * q_dot;
    M_inv = inv(M);

    tau = pinv(J * M_inv) * (x_ddot + J * M_inv * n);
    %tau = pinv(J * M_inv) * x_ddot + n;
    if nargin > 6
        J_inv = pinv(J);
        num_rows = size(J_inv, 1);
        I = eye(num_rows);
        H = (inv(J * M_inv * J') * J * M_inv); %generalized inverse
        tau = tau + M * (I - J' * H) * tau_0;
    end

    tau = simplify(tau);
end

%τ₀ = (I - (JM⁻¹)⁺(JM⁻¹))v = M(I - J⁺J)v 