function Y = extract_Y(Ya, dynamic_symbols)
    % Crediti a Gabbo (thx <3)
    % M: inertia matrix
    % C: christoffel matrix
    % G: potential enrgy matrix
    % joint_acceleration: array with joint acceleration
    % this four items are use to build dynamic model

    % Ya: in the form M*ddq + C + G (eventually + Fv*dq) with a1.. a_n already included
    % dynamic_symbols: how substitute coefficient, so [a1 a2 ... a_n]

    Y = sym([]);
    Y_columns = length(dynamic_symbols);
    dynamic_symbols = reshape(dynamic_symbols, 1, Y_columns);

    for i=1:Y_columns
        
        subs_array = zeros(1, Y_columns);
        subs_array(i) = 1;

        Yi = subs(Ya, dynamic_symbols, subs_array);

        Y = [Y, Yi];
    end

    Y = simplify(Y);

    %check if Y * a = Ya
    %assert(Y * dynamic_symbols == Ya)

    
end

