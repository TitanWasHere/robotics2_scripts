function vc = computeCoMVelocity(w, v, r_com)
% COMPUTECOMVELOCITY Calcola la velocità del centro di massa (vc) dati w, v e r_com.
% Input:
%   w - Cell array delle velocità angolari (i_wi) per ogni frame
%   v - Cell array delle velocità lineari (i_vi) per ogni frame
%   r_com - Cell array delle posizioni del CoM (i_r_ci) per ogni frame
% Output:
%   vc - Cell array delle velocità CoM (i_vci) per ogni frame

n = length(w);
vc = cell(1, n);

for i = 1:n
    % Velocità del CoM: vc = v + w × r_com
    vc{i} = simplify(v{i} + cross(w{i}, r_com{i}));
end
end