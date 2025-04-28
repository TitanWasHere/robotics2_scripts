function [CoM_mobile] = calculate_CoM_mobile(DH_table, dc)
% CALCOLATE_COM_MOBILE_SYMBOLIC Calcola la posizione simbolica del CoM
% Input:
%   DH_table - Tabella parametri DH simbolica [alpha a d theta]
%   dc - Vettore delle distanze simboliche dei CoM [dc1; dc2; ...; dcn]
% Output:
%   CoM_mobile - Cell array con i vettori posizione CoM simbolici

    n_links = size(DH_table,1);
    CoM_mobile = cell(1,n_links);
    
    % Estrai parametri DH (simbolici)
    a = DH_table(:,2);  % Vettore dei parametri a
    
    for i = 1:n_links
        if i == 1
            % CoM Link 1 (speciale)
            CoM_mobile{i} = [0; dc(i); 0];
        else
            % CoM Link i>1 (formula generale)
            CoM_mobile{i} = [dc(i) - a(i); 0; 0];
        end
    end
    
    % Visualizzazione risultati simbolici
    fprintf('Posizioni CoM simboliche nei frame mobili:\n');
    for i = 1:n_links
        fprintf('Link %d:\n', i);
        disp(CoM_mobile{i})
    end
end