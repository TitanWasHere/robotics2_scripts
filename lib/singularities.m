function[Js] = singularities(J, joints)
    [n_rows, n_col] = size(J);
    if n_rows ~= n_col
        fprintf("[WARNING]: matrix not squared, calulating J(q) * J(q)^T\n")
        J_new = J * transpose(J);
        det_J = simplify(det(J_new));
        fprintf("[det]: %s\n", det_J);
        if det_J == 0
            J_new = transpose(J) * J; % if 0 then transpose(J)*J
            det_J = simplify(det(J_new));
            fprintf("[det]: %s\n", det_J);
        end
    else
        det_J = simplify(det(J));
        fprintf("[det]: %s\n", det_J);
    end

    
    fprintf("NOW YOU SHOULD MODIFY IN ORDER TO SUBSTITUTE TO THE JOINT ITS SINGULARITY")
    pause
    % es.
    Js = subs(J, [joints(2)], [pi/2]);
    display(Js)
    
    % dim. of null space
    dim_null = n_col - rank(Js);
    fprintf("[n]: %d\n", n_col);
    fprintf("[rank]: %d\n", rank(Js))
    fprintf("[dim_null]: rank-n = %d\n", dim_null)

    % null space --> vel. ai joint che danno vel 0 all'ee
    null_space = simplify(null(Js));
    fprintf("[NULL SPACE]:")
    display(null_space);

    % dim. of range space:
    dim_range = rank(Js);
    display(dim_range);
    fprintf("[RANGE SPACE]:")
    display(orth(Js))
    
    % TO FIND RANGE SPACE: TAKE THE 'dim_range' COLS MORE CONVENIENT AND WE
    % PUT THEM INTO THE RANGE SPACE (maybe with R(J) = span { (...), (...)}
    % --> vel. che vengono mandate dai joint all'ee
end