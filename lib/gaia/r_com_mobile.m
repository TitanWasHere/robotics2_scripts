function [r_com] = r_com_mobile(DH_table)
% CALCOLATE_COM_MOBILE Calcola la posizione del CoM nei frame mobili
% Input:
%   DH_table - Tabella parametri DH [alpha a d theta]
% Output:
%   r_com - Cell array con i vettori posizione CoM nei frame mobili

    n_links = size(DH_table,1);
    r_com = cell(1,n_links);
    
    % Crea vettore simbolico delle distanze CoM
    dc = sym('dc', [n_links 1], 'real');  % [dc1; dc2; ...; dcn]
    
    % Estrai parametri a dalla tabella DH
    a = DH_table(:,2);  % Vettore dei parametri a
    
    for i = 1:n_links
        if i == 1
            % CoM Link 1 (speciale)
            r_com{i} = [0; dc(i); 0];
        else
            % CoM Link i>1 (formula generale)
            r_com{i} = [dc(i) - a(i); 0; 0];
        end
    end
    
    % Visualizzazione risultati (opzionale)
    fprintf('Posizioni CoM nei frame mobili:\n');
    for i = 1:n_links
        fprintf('Link %d:\n', i);
        disp(r_com{i})
    end
end