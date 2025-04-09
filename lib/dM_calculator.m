function dM = dM_calculator(M, q, dq)
    % Input:
    %   - M: inertia matrix
    %   - q: vertical vector of joints
    %   - dq: vertical vector of joint velocities
    % NOTE:
    %   Joints must be ordered by ascending order: (q1, q2, ... )
    %
    % Output:
    %   - Time derivative of M

    n = size(q,1);
    dM = zeros(size(M,1), size(M,2));
    for i = 1:n
        dM = dM + diff(M, q(i))*dq(i);
    end
end