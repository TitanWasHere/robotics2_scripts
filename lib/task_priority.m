function q_dot = task_priority(J1, r1_dot, J2, r2_dot, v2)

    J1_inv = pinv(J1);
    J2_inv = pinv(J2);

    P1 = proj_NullSpace(J1);
    P2 = proj_NullSpace(J2);

    if nargin < 5
        %v2 is the same legth as the rows of J1_inv
        v2 = zeros(size(J1_inv, 1), 1);
    end

    J2P1 = simplify(J2 * P1);
    J1_inv_r1_dot = simplify(J1_inv * r1_dot);
    J2P1_inv = pinv(J2P1);

    %q_dot = J1_inv_r1_dot + J2P1_inv * (r2_dot - J2  * J1_inv_r1_dot) + P1 * (eye(size(J2P1, 2)) - J2P1_inv * J2P1) * v2;
    q_dot = J1_inv_r1_dot + P1 * J2P1_inv * (r2_dot - J2  * J1_inv_r1_dot) + P1 * (eye(size(J2P1, 2)) - J2P1_inv * J2P1) * v2;
    %q_dot = J1_inv_r1_dot + P1 * J2P1_inv * (r2_dot - J2  * J1_inv_r1_dot) + P1 * proj_NullSpace(J2P1) * v2;

    q_dot = simplify(q_dot);