function [x, iter] = Jacobi(A, b, x_initial, maxiter, tol)
    x = x_initial;
    Dinv = 1 ./ diag(A); % Vector of diagonal inverses
    normb = norm(b, 2);
    for iter = 1:maxiter
        r = b - A * x;
        x_new = x + Dinv .* r;
        if norm(b - A * x_new, 2) <= tol * normb
            x = x_new;
            return;
        end
        x = x_new;
    end
end