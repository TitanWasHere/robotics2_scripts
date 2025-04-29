function res = proj_NullSpace(J, W)

    if nargin == 2
        I = eye(length(J));
        res = (I - weighted_pinv(J, W)*J);
        res = simplify(res);
        return;

    end

    I = eye(length(J));
    res = (I-pinv(J)*J);
    res = simplify(res);
end