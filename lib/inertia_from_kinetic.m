function M = inertia_from_kinetic(T, qdot)
    % Takes as inputs:
    %   - T = the kinetic energy of the robot 
    %   - qdot_vector = the vector of q_dot ex: [qd1,qd2,qd3]
    %
    % Output:
    %   - M = the inertia matrix 

    M = hessian(T, qdot);
end