function [x, iter] = SOR(omega, A, b, x_initial, maxiter, tol)
    x = x_initial;
    n = length(b);
    normb = norm(b, 2);
    A_diag = diag(A); % Precompute diagonal

    [row, col, val] = find(A); % Get sparse structure

    rowcols = accumarray(row, col, [n, 1], @(c) {c});
    rowvals = accumarray(row, val, [n, 1], @(v) {v});

    for iter = 1:maxiter
        x_old = x;
        for i = 1:n
            cols = rowcols{i};
            vals = rowvals{i};
            % Exclude diagonal
            mask_left = cols < i;
            mask_right = cols > i;
            sum1 = vals(mask_left)' * x(cols(mask_left));
            sum2 = vals(mask_right)' * x_old(cols(mask_right));
            x(i) = (1 - omega) * x_old(i) + omega / A_diag(i) * (b(i) - sum1 - sum2);
        end
        r = b - A * x;
        if norm(r, 2) <= tol * normb
            return;
        end
    end
end