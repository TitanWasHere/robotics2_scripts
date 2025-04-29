function dq = projected_gradient(J, q, v, H)
    % J: Jacobian of the end effector position vector
    % q: joints vector [q1,q2,...]
    % v: end effector velocity vector [vx; vy; vz];
    % H:
    n = length(J);
    I = eye(n);
    dq = pinv(J) * v + (I - pinv(J) * J) * gradient(H, q);



end