function [x, iter] = GS(A, b, x_initial, maxiter, tol)
    x = x_initial;
    normb = norm(b, 2);
    L = tril(A);
    U = triu(A,1);
    for iter = 1:maxiter
        x_new = L \ (b - U * x);
        if norm(b - A * x_new, 2) <= tol * normb
            x = x_new;
            return;
        end
        x = x_new;
    end
end