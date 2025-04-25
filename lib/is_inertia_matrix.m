function is_an_inertia_matrix = is_inertia_matrix(matrix)
%IS_INERTIA_MATRIX Validate if a matrix can be an inertia matrix.
%   Checks symmetry and positive definiteness, with verbose output.

    is_an_inertia_matrix = false;  % Default

    % 1) Symmetry Check
    if isequal(matrix, matrix.')
        fprintf('✅ Matrix is symmetric.\n');
    else
        fprintf('❌ Matrix is NOT symmetric.\n');
        return;
    end

    % 2) Determinant Check
    if isa(matrix, 'sym')
        fprintf('⚠ Matrix is symbolic. Please verify the determinant manually.\n');
        fprintf('Determinant: ');
        disp(simplify(det(matrix)));
    else
        if det(matrix) > 0
            fprintf('✅ Determinant is positive.\n');
        else
            fprintf('❌ Determinant is non-positive.\n');
            fprintf('🔴 Invalid inertia matrix (determinant ≤ 0).\n');
            return;
        end
    end

    % 3) Eigenvalue Check
    eigenvalues = eig(matrix);
    if isa(matrix, 'sym')
        fprintf('⚠ Matrix is symbolic. Please verify eigenvalues manually.\n');
        fprintf('   eigenvalues: \n');
        disp(simplify(eigenvalues));
        is_an_inertia_matrix = [];
    elseif all(eigenvalues > 0)
        fprintf('✅ All eigenvalues are positive.\n');
        is_an_inertia_matrix = true;
    else
        fprintf('❌ Some eigenvalues are non-positive.\n');
        fprintf('🔴 Invalid inertia matrix (not positive definite).\n');
        return;
    end

    % Final suggestion
    if isequal(is_an_inertia_matrix, true)
        fprintf('\n🟢 Valid inertia matrix.\n');
        fprintf('🔍 Reminder: Make sure the matrix does not depend directly on q1 (arbitrary coordinate).\n');
    end
end
