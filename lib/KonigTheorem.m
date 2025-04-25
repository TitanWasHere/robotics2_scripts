function T = KonigTheorem(m, vc, w, I)
    % Takes as inputs:
    %   - m = the mass of the robot
    %   - vc = the velocity of the center of mass
    %   - I = the inertia tensor
    %
    % Output:
    %   - T = the kinetic energy of the robot 
    %%
    %   T = KonigTheorem(m, vc) computes translational kinetic energy only.
    %   T = KonigTheorem(m, vc, w, I) computes full kinetic energy including rotation.
    %%
    T = 0.5 * m * simplify(norm(vc)^2);

    % Check if rotational terms are provided
    if nargin == 4
        T = T + 0.5 * w' * I * w;
    elseif nargin ~= 2
        error('❌ Function expects either 2 or 4 input arguments.');
    end
end 