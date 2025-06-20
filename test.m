m = 16;
alpha = 6.4e-2;
K = 3;
tol = 1e-2;
maxiter = 1000;
methods = {@Jacobi, @GS, @CG};
method_names = {'Jacobi', 'Gauss-Seidel', 'Conjugate Gradient'};
meth = 1; % Choose 1 for Jacobi, 2 for GS, 3 for CG

X = set_image(m);
u0 = FormRHS(X);
n = m^2;
x0 = zeros(n,1);

A = FormMatrix(x0, alpha); % Only once, outside the K loop
x = x0;
total_iter = 0;
for k = 1:K
    b = u0;
    [x, iter] = methods{meth}(A, b, x, maxiter, tol);
    total_iter = total_iter + iter;
end

% Show result
figure; imagesc(reshape(x, m, m)); colormap gray; axis image; title(method_names{meth});
disp(['Total iterations: ', num2str(total_iter)]);