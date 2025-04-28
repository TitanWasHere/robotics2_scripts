function [A_consec, R, r_pre, r_mov] = my_DH_parser(DH_table)
% Funzione per calcolare le matrici omogenee consecutive (A_consec),
% le matrici di rotazione (R), i vettori posizione nei frame precedenti (r_pre),
% e nei frame mobili (r_mov) a partire dalla tabella dei parametri DH.
%
% Input:
%   DH_table - Tabella dei parametri DH [alpha, a, d, theta]
% Output:
%   A_consec - Matrici omogenee consecutive (frame i-1 → frame i)
%   R        - Matrici di rotazione (frame i-1 → frame i)
%   r_pre    - Posizione origine frame i nel frame precedente (i-1)
%   r_mov    - Posizione origine frame i-1 nel frame mobile (i)

rows = size(DH_table, 1);
A_consec = cell(1, rows);
R = cell(1, rows);
r_pre = cell(1, rows);  % Nome modificato da r_consec
r_mov = cell(1, rows);  % Nome modificato da r_moving

for i = 1:rows
    alpha = DH_table(i,1);
    a     = DH_table(i,2);
    d     = DH_table(i,3);
    theta = DH_table(i,4);

    % Matrice omogenea i-1 → i
    A = [cos(theta), -cos(alpha)*sin(theta),  sin(alpha)*sin(theta), a*cos(theta);
         sin(theta),  cos(alpha)*cos(theta), -sin(alpha)*cos(theta), a*sin(theta);
         0,           sin(alpha),             cos(alpha),            d;
         0,           0,                      0,                     1];

    A = simplify(A);
    A_consec{i} = A;
    
    % Estrazione componenti
    R{i} = simplify(A(1:3,1:3));       % Rotazione i-1 → i
    r_pre{i} = simplify(A(1:3,4));     % Posizione O_i nel frame i-1 (precedente)
    
    A_inv = simplify(inv(A));
    r_mov{i} = simplify(-A_inv(1:3,4)); % Posizione O_{i-1} nel frame i (mobile)
end
end