function T = KonigTheorem(m, vc, w, I)
    % Takes as inputs:
    %   - m = the mass of the robot
    %   - vc = the velocity of the center of mass
    %   - I = the inertia tensor
    %
    % Output:
    %   - T = the kinetic energy of the robot 

    T = 0.5 * m * simplify(norm(vc)^2) + 0.5 * w' * I * w;
end 