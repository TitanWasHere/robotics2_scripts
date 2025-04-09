function S = S_calculator(C, dq)
    % Function that calculate the S matrix given the CHRISTOFFEL Symobls
    % Input
    %   - C: C{i} is the i-th Christoffel matrix
    %   - dq: joint velocities vertical vector
    %
    % Output:
    %   - S: matrix S

    n = size(dq,1);
    S = [];
    for i=1:n
        Si = dq' * C{i};
        S = [S; Si];
    end
end