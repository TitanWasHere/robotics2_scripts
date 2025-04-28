function I = InertiaFromCom(r_com)
% InertiaFromCom - Genera matrici di inerzia simboliche per ogni link
%
% Syntax: I = InertiaFromCom(r_com)
%
% Inputs:
%    r_com - cell array {r1, r2, ..., rn} con i vettori posizione centro di massa [x, y, z]
%
% Outputs:
%    I - cell array delle matrici di inerzia simboliche {I1, I2, ..., In}

    % Inizializza cell array output
    I = cell(1, length(r_com));

    for i = 1:length(r_com)
        % Definizione dei simboli per il link i
        Ic_ix = sym(sprintf('Ic_%dx', i), 'real');
        Ic_iy = sym(sprintf('Ic_%dy', i), 'real');
        Ic_iz = sym(sprintf('Ic_%dz', i), 'real');

        % Estrai posizione centro di massa
        r = r_com{i};

        % Se r è simbolico, assumiamo che solo una componente sia ≠ 0
        % Allora controlliamo direttamente:
        is_x = ~isequal(r(1), sym(0)) && isequal(r(2), sym(0)) && isequal(r(3), sym(0));
        is_y = ~isequal(r(2), sym(0)) && isequal(r(1), sym(0)) && isequal(r(3), sym(0));
        is_z = ~isequal(r(3), sym(0)) && isequal(r(1), sym(0)) && isequal(r(2), sym(0));

        % Gestione simmetria
        if is_x
            % Asse principale X
            I{i} = diag([Ic_ix, Ic_iz, Ic_iz]);
        elseif is_y
            % Asse principale Y
            I{i} = diag([Ic_ix, Ic_iy, Ic_ix]);
        elseif is_z
            % Asse principale Z
            I{i} = diag([Ic_ix, Ic_ix, Ic_iz]);
        else
            % Nessuna simmetria speciale
            I{i} = diag([Ic_ix, Ic_iy, Ic_iz]);
        end
    end

end
