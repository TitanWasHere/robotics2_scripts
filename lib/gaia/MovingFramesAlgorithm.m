function [w, v, vc] = MovingFramesAlgorithm(R, r_pre, r_com, q_dot, sigma)
% MOVINGFRAMESALGORITHM Calcola w, v e vc per tutti i frame
% Input:
%   R - Cell array di matrici di rotazione (i-1_R_i)
%   r_pre - Cell array di vettori posizione nel frame precedente (i-1_r_i)
%   r_com - Cell array di vettori posizione CoM nel frame mobile (i_r_ci)
%   q_dot - Vettore delle velocità dei giunti [q_dot1; q_dot2; ...]
%   sigma - Vettore tipo giunto (0: revolute, 1: prismatico)
% Output:
%   w - Cell array delle velocità angolari (i_wi)
%   v - Cell array delle velocità lineari (i_vi)
%   vc - Cell array delle velocità CoM (i_vci)

n = length(R);
w = cell(1, n);
v = cell(1, n);
vc = cell(1, n);

% Inizializzazione: Frame 0 (base)
w_prev = [0; 0; 0];
v_prev = [0; 0; 0];

for i = 1:n
    R_prev_to_current = R{i}.'; % i-1_R_i trasposto = i_R_i-1
    z_prev = [0; 0; 1]; % Asse z in RF_{i-1}

    % 1. Calcolo velocità angolare i_wi
    w_current = R_prev_to_current * (w_prev + (1 - sigma(i)) * q_dot(i) * z_prev);
    w{i} = simplify(w_current);

    % 2. Calcolo velocità lineare i_vi
    cross_w_r_pre = cross(w{i}, r_pre{i});
    v_current = R_prev_to_current * (v_prev + sigma(i) * q_dot(i) * z_prev + cross_w_r_pre);
    v{i} = simplify(v_current);

    % 3. Calcolo velocità del centro di massa i_vci
    cross_w_r_com = cross(w_current, r_com{i});
    vc_current = v_current + cross_w_r_com;
    vc{i} = simplify(vc_current);

    % Aggiorna valori per il prossimo link
    w_prev = w_current;
    v_prev = v_current;
end
end
