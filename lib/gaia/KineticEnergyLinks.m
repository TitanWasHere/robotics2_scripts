function [T, Ttot] = KineticEnergyLinks(m, vc , w, I)
    % KineticEnergyLinks - Calcola energia cinetica per ogni link e totale
    %
    % Syntax: [T, Ttot] = KineticEnergyLinks(m, v, w, I)
    %
    % Inputs:
    %    m - vettore delle masse [m1 m2 ... mn]
    %    v - cell array delle velocità lineari {v1, v2, ..., vn}
    %    w - cell array delle velocità angolari {w1, w2, ..., wn}
    %    I - cell array delle matrici di inerzia {I1, I2, ..., In}
    %
    % Outputs:
    %    T - cell array dell'energia cinetica di ciascun link {T1, T2, ..., Tn}
    %    Ttot - energia cinetica totale del sistema

    % Verifica consistenza input
    n = length(m);
    if length(vc) ~= n || length(w) ~= n || length(I) ~= n
        error('Tutti gli input devono avere la stessa lunghezza');
    end
    
    T = cell(1, n); % Inizializza cell array per le energie
    
    % Calcola energia cinetica per ogni link
    for i = 1:n
        % Verifica dimensioni
        if numel(vc{i}) ~= 3 || numel(w{i}) ~= 3 || ~isequal(size(I{i}), [3 3])
            error('Input non validi: v e w devono essere 3x1, I deve essere 3x3');
        end
        
        % Componente traslazionale
        T_trans = (1/2)*m(i)*(vc{i}'*vc{i});
        
        % Componente rotazionale
        T_rot = (1/2)*w{i}'*I{i}*w{i};
        
        % Energia totale del link
        T{i} = simplify(T_trans + T_rot);
    end
    
    % Energia cinetica totale
    Ttot = sym(0);
    for i = 1:n
        Ttot = Ttot + T{i};
    end
    Ttot = simplify(Ttot);
end