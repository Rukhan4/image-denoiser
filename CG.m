function [x, iter] = CG(A, b, x_initial, maxiter, tol)
    x = x_initial;
    r = b - A * x;
    p = r;
    normb = norm(b, 2);

    for iter = 1:maxiter
        Ap = A * p;
        alpha = (r' * r) / (p' * Ap);
        x = x + alpha * p;
        r_new = r - alpha * Ap;

        if norm(r_new, 2) <= tol * normb
            return;
        end

        beta = (r_new' * r_new) / (r' * r);
        p = r_new + beta * p;
        r = r_new;
    end
end
