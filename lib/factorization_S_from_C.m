function S_s = factorization_S_from_C(cell_C_k, q, q_dot, print_info)
    %Inputs:
    %   cell_C_k : cell array of Ck
    %   q : (Nx1) vector of joint configuration
    %   q_dot : (Nx1) vector of joint velocity
    %   print_info : 'true' to print info
    %Output:
    %   S : skew-sym matrix

    S_s = sym([]);
    
    n = length(q); 
    for i=1:n
        s_i = q_dot'*cell_C_k{i};
        S_s = [S_s; s_i];
    end

    if nargin < 4
        print_info = false;
    end
    
    if print_info
        fprintf("Row i-th of S matrix is given by s_i = q_dot' * c_i\n");
        fprintf("Result S:\n");
        display(S_s);
        fprintf("----------------------------------------------------\n");
    end
end