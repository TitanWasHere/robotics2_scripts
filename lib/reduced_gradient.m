function [q_dot, reduced_gradient] = reduced_gradient(Ja, Jb, H, v, T, q)
%REDUCEDGRADIENT  compute reduced gradient and joint velocities
%
% [q_dot, g_red] = reducedGradient(Ja, Jb, H, r_dot, T, q)
%
% Inputs:
%   Ja    – m×m active‐Jacobian
%   Jb    – m×(n–m) passive‐Jacobian
%   H     – symbolic scalar Hamiltonian, function of q
%   r_dot – m×1 task‐space velocity
%   T     – (n–m)×(n–m) “identity” matrix (usually eye(n–m))
%   q     – 1×n symbolic joint vector [q1,q2,…,qn]
%
% Outputs:
%   q_dot – n×1 joint‐velocity vector [q̇a; q̇b]
%   g_red – (n–m)×1 reduced gradient 
%
% Example:
%   syms q1 q2 q3 real
%   q  = [q1 q2 q3];
%   Ja = sym(eye(2));            % e.g. m=2
%   Jb = sym([1;2]);             % 2×1
%   H  = q1^2+sin(q2)+q3;        % scalar Hamiltonian
%   r_dot = sym([.1; .2]);
%   T = eye(1);
%   [q_dot, g_red] = reducedGradient(Ja,Jb,H,r_dot,T,q);
    invJa = inv(Ja);
    %q_dot = [invJa ; 0] * v + [ - invJa*Jb ; T ] * [-(invJa * Jb)' , T ] * gradient(H, q);
    grad_H = gradient(H, q);
    reduced_gradient = [ -(invJa * Jb)' , 1] * T * grad_H;
    q_dot = T' * [invJa * (v - Jb* reduced_gradient) ; reduced_gradient];
end
