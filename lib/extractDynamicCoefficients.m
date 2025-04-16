function [dynamic_params, M_regrouped] = extractDynamicCoefficients(M, q)
% EXTRACTDYNAMICCOEFFICIENTS Estrae i parametri dinamici dalla matrice di inerzia
% e riorganizza la matrice in una forma parametrizzata
%
% Input:
%   M: Matrice di inerzia simbolica
%   q: Vettore delle variabili di giunto [q1, q2, ..., qn]
%
% Output:
%   dynamic_params: Struttura contenente i parametri dinamici identificati
%   M_regrouped: Matrice di inerzia riorganizzata in termini dei parametri dinamici
%
% Esempio di utilizzo:
%   syms q1 q2 real
%   q = [q1, q2];
%   % Dopo aver calcolato M...
%   [params, M_param] = extractDynamicCoefficients(M, q);

    % Verifica degli input
    if nargin < 2
        error('Sono richiesti due argomenti: matrice di inerzia e variabili di giunto');
    end
    
    % Dimensione della matrice di inerzia
    [n, ~] = size(M);
    
    % Verifica che M sia quadrata
    if n ~= size(M, 2)
        error('La matrice di inerzia deve essere quadrata');
    end
    
    % Inizializza le strutture per i parametri dinamici
    dynamic_params = struct('values', struct());
    param_symbols = struct();
    
    % Crea simboli per i termini trigonometrici
    for i = 1:n
        s_name = sprintf('s%d', i);
        c_name = sprintf('c%d', i);
        dynamic_params.(s_name) = sin(q(i));
        dynamic_params.(c_name) = cos(q(i));
    end
    
    % Inizializza la matrice riorganizzata
    M_regrouped = sym(zeros(n, n));
    
    % Contatore per i coefficienti
    coef_count = 1;
    
    % Analizza ogni elemento della matrice
    for i = 1:n
        for j = i:n  % Sfrutta la simmetria della matrice di inerzia
            element = M(i, j);
            
            % Salta se l'elemento è zero
            if element == 0
                continue;
            end
            
            % Cerca dipendenze da termini trigonometrici
            has_trig = false;
            for k = 1:n
                if has(element, sin(q(k))) || has(element, cos(q(k)))
                    has_trig = true;
                    break;
                end
            end
            
            % Se non ha termini trigonometrici, è un elemento costante
            if ~has_trig
                param_name = sprintf('a%d', coef_count);
                param_symbols.(param_name) = sym(param_name, 'real');
                dynamic_params.values.(param_name) = element;
                M_regrouped(i, j) = param_symbols.(param_name);
                if i ~= j
                    M_regrouped(j, i) = M_regrouped(i, j);  % Simmetria
                end
                coef_count = coef_count + 1;
            else
                % Analisi per termini trigonometrici
                terms = struct();
                
                % Termine costante
                terms.constant = 0;
                
                % Termini con sin
                for k = 1:n
                    s_name = sprintf('s%d', k);
                    terms.(s_name) = 0;
                end
                
                % Termini con cos
                for k = 1:n
                    c_name = sprintf('c%d', k);
                    terms.(c_name) = 0;
                end
                
                % Termini prodotto sin*cos
                for k = 1:n
                    for l = 1:n
                        sc_name = sprintf('s%dc%d', k, l);
                        terms.(sc_name) = 0;
                    end
                end
                
                % Termini prodotto sin*sin
                for k = 1:n
                    for l = k:n
                        ss_name = sprintf('s%ds%d', k, l);
                        terms.(ss_name) = 0;
                    end
                end
                
                % Termini prodotto cos*cos
                for k = 1:n
                    for l = k:n
                        cc_name = sprintf('c%dc%d', k, l);
                        terms.(cc_name) = 0;
                    end
                end
                
                % Estrai i termini costanti
                try
                    % Se ci sono variabili, cerca di estrarre la parte costante
                    vars = symvar(element);
                    if ~isempty(vars)
                        [c, t] = coeffs(element, vars);
                        for k = 1:length(t)
                            if t(k) == 1  % Termine costante
                                terms.constant = terms.constant + c(k);
                            end
                        end
                    else
                        % Se non ci sono variabili, l'intero elemento è costante
                        terms.constant = element;
                    end
                catch
                    % Ignora errori nell'estrazione
                end
                
                % Estrai i termini con seno e coseno singoli
                for k = 1:n
                    % Estrai coefficiente di sin(q_k)
                    try
                        s_term = sin(q(k));
                        if has(element, s_term)
                            s_coef = coeffs_for_term(element, s_term);
                            if ~isempty(s_coef) && s_coef ~= 0
                                s_name = sprintf('s%d', k);
                                terms.(s_name) = s_coef;
                            end
                        end
                    catch
                        % Ignora errori
                    end
                    
                    % Estrai coefficiente di cos(q_k)
                    try
                        c_term = cos(q(k));
                        if has(element, c_term)
                            c_coef = coeffs_for_term(element, c_term);
                            if ~isempty(c_coef) && c_coef ~= 0
                                c_name = sprintf('c%d', k);
                                terms.(c_name) = c_coef;
                            end
                        end
                    catch
                        % Ignora errori
                    end
                end
                
                % Estrai prodotti sin*cos
                for k = 1:n
                    for l = 1:n
                        try
                            sc_term = sin(q(k))*cos(q(l));
                            if has(element, sc_term)
                                sc_coef = coeffs_for_term(element, sc_term);
                                if ~isempty(sc_coef) && sc_coef ~= 0
                                    sc_name = sprintf('s%dc%d', k, l);
                                    terms.(sc_name) = sc_coef;
                                end
                            end
                        catch
                            % Ignora errori
                        end
                    end
                end
                
                % Estrai prodotti sin*sin
                for k = 1:n
                    for l = k:n
                        try
                            ss_term = sin(q(k))*sin(q(l));
                            if has(element, ss_term)
                                ss_coef = coeffs_for_term(element, ss_term);
                                if ~isempty(ss_coef) && ss_coef ~= 0
                                    ss_name = sprintf('s%ds%d', k, l);
                                    terms.(ss_name) = ss_coef;
                                end
                            end
                        catch
                            % Ignora errori
                        end
                    end
                end
                
                % Estrai prodotti cos*cos
                for k = 1:n
                    for l = k:n
                        try
                            cc_term = cos(q(k))*cos(q(l));
                            if has(element, cc_term)
                                cc_coef = coeffs_for_term(element, cc_term);
                                if ~isempty(cc_coef) && cc_coef ~= 0
                                    cc_name = sprintf('c%dc%d', k, l);
                                    terms.(cc_name) = cc_coef;
                                end
                            end
                        catch
                            % Ignora errori
                        end
                    end
                end
                
                % Assegna parametri ai termini trovati
                expr = 0;
                
                % Termine costante
                if terms.constant ~= 0
                    param_name = sprintf('a%d', coef_count);
                    param_symbols.(param_name) = sym(param_name, 'real');
                    dynamic_params.values.(param_name) = terms.constant;
                    expr = expr + param_symbols.(param_name);
                    coef_count = coef_count + 1;
                end
                
                % Termini con sin
                for k = 1:n
                    s_name = sprintf('s%d', k);
                    if isfield(terms, s_name) && terms.(s_name) ~= 0
                        param_name = sprintf('a%d', coef_count);
                        param_symbols.(param_name) = sym(param_name, 'real');
                        dynamic_params.values.(param_name) = terms.(s_name);
                        expr = expr + param_symbols.(param_name) * dynamic_params.(s_name);
                        coef_count = coef_count + 1;
                    end
                end
                
                % Termini con cos
                for k = 1:n
                    c_name = sprintf('c%d', k);
                    if isfield(terms, c_name) && terms.(c_name) ~= 0
                        param_name = sprintf('a%d', coef_count);
                        param_symbols.(param_name) = sym(param_name, 'real');
                        dynamic_params.values.(param_name) = terms.(c_name);
                        expr = expr + param_symbols.(param_name) * dynamic_params.(c_name);
                        coef_count = coef_count + 1;
                    end
                end
                
                % Termini prodotto sin*cos
                for k = 1:n
                    for l = 1:n
                        sc_name = sprintf('s%dc%d', k, l);
                        if isfield(terms, sc_name) && terms.(sc_name) ~= 0
                            param_name = sprintf('a%d', coef_count);
                            param_symbols.(param_name) = sym(param_name, 'real');
                            dynamic_params.values.(param_name) = terms.(sc_name);
                            expr = expr + param_symbols.(param_name) * dynamic_params.(sprintf('s%d', k)) * dynamic_params.(sprintf('c%d', l));
                            coef_count = coef_count + 1;
                        end
                    end
                end
                
                % Termini prodotto sin*sin
                for k = 1:n
                    for l = k:n
                        ss_name = sprintf('s%ds%d', k, l);
                        if isfield(terms, ss_name) && terms.(ss_name) ~= 0
                            param_name = sprintf('a%d', coef_count);
                            param_symbols.(param_name) = sym(param_name, 'real');
                            dynamic_params.values.(param_name) = terms.(ss_name);
                            expr = expr + param_symbols.(param_name) * dynamic_params.(sprintf('s%d', k)) * dynamic_params.(sprintf('s%d', l));
                            coef_count = coef_count + 1;
                        end
                    end
                end
                
                % Termini prodotto cos*cos
                for k = 1:n
                    for l = k:n
                        cc_name = sprintf('c%dc%d', k, l);
                        if isfield(terms, cc_name) && terms.(cc_name) ~= 0
                            param_name = sprintf('a%d', coef_count);
                            param_symbols.(param_name) = sym(param_name, 'real');
                            dynamic_params.values.(param_name) = terms.(cc_name);
                            expr = expr + param_symbols.(param_name) * dynamic_params.(sprintf('c%d', k)) * dynamic_params.(sprintf('c%d', l));
                            coef_count = coef_count + 1;
                        end
                    end
                end
                
                % Assegna l'espressione completa alla matrice riorganizzata
                M_regrouped(i, j) = expr;
                if i ~= j
                    M_regrouped(j, i) = M_regrouped(i, j);  % Simmetria
                end
            end
        end
    end
    
    % Aggiungi tutti i simboli dei parametri alla struttura di output
    fields = fieldnames(param_symbols);
    for i = 1:length(fields)
        dynamic_params.(fields{i}) = param_symbols.(fields{i});
    end
    
    % Stampa l'interpretazione fisica dei parametri dinamici
    disp('Interpretazione fisica dei parametri dinamici:');
    fprintf('Per un manipolatore con %d giunti:\n', n);
    
    fields = fieldnames(param_symbols);
    for i = 1:length(fields)
        fprintf('%s: Parametro dinamico che può rappresentare massa, inerzia o combinazioni\n', fields{i});
    end
    
    % Suggerimenti su come interpretare i parametri
    disp('');
    disp('Interpretazione tipica dei parametri per un manipolatore a 2 giunti:');
    disp('a1: Combinazione di inerzie dei link e masse (I₁ + m₁r₁² + I₂ + m₂l₁² + m₂r₂²)');
    disp('a2: Accoppiamento inerziale tra i giunti (m₂l₁r₂)');
    disp('a3: Coefficiente che moltiplica sin(q₂) (spesso legato a effetti gravitazionali)');
    disp('a4: Inerzia del secondo link (I₂ + m₂r₂²)');
    
    % Esempio di interpretazione fisica della matrice risultante
    disp('');
    disp('Per un manipolatore a 2 giunti, la matrice di inerzia tipicamente ha forma:');
    disp('M(q) = [ a₁ + 2a₂c₂, a₄ + a₂c₂ ]');
    disp('       [ a₄ + a₂c₂,       a₄    ]');
    disp('dove c₂ = cos(q₂)');
end

function coef = coeffs_for_term(expr, term)
    % Estrae il coefficiente per un termine specifico
    try
        % Cerca di usare coeffs direttamente
        [c, t] = coeffs(expr, term);
        idx = find(t == term);
        if ~isempty(idx)
            coef = c(idx);
        else
            % Prova con la derivata parziale
            coef = diff(expr, term);
            % Se la derivata contiene ancora il termine, non è lineare
            if has(coef, term)
                coef = [];
            end
        end
    catch
        % In caso di errore, prova con l'approccio della sostituzione
        try
            % Sostituisci il termine con una variabile simbolica ausiliaria
            aux_var = sym('aux_var', 'real');
            expr_subst = subs(expr, term, aux_var);
            
            % Estrai il coefficiente lineare
            [c, t] = coeffs(expr_subst, aux_var);
            idx = find(t == aux_var);
            if ~isempty(idx)
                coef = c(idx);
            else
                coef = [];
            end
        catch
            % Se tutto fallisce, restituisci vuoto
            coef = [];
        end
    end
end