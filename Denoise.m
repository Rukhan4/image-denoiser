clear; clc;
profile on;

% Parameters
sizes = [16, 32, 64]; 
alphas = [6.4e-2, 3.2e-2, 1.6e-2];
K = 8; 
tol = 1e-2;
maxiter = 20000;

methods = {@Jacobi, @GS, @CG}; 
method_names = {'Jacobi', 'Gauss-Seidel', 'Conjugate Gradient'};
omega_opt = [1.12, 1.10, 1.08]; 

results_time = zeros(length(sizes), length(methods)+1);
results_iter = zeros(length(sizes), length(methods)+1);

for s = 1:length(sizes)
    m = sizes(s);
    n = m^2;
    alpha = alphas(s);
   
    X = set_image(m); 
    u0 = FormRHS(X);

   
    x0 = zeros(n,1);

    % Loop over methods
    for meth = 1:length(methods)
        x = x0;
        total_iter = 0;
        tic;
        for k = 1:K
            disp("Running method: " + method_names{meth} + " for size: " + m + " run: " + k)
            A = FormMatrix(x, alpha);
            b = u0;
            [x, iter] = methods{meth}(A, b, x, maxiter, tol);
            total_iter = total_iter + iter;
        end
        elapsed = toc;
        results_time(s, meth) = elapsed;
        results_iter(s, meth) = total_iter;
        if s == 1
           figure; imagesc(reshape(x, m, m)); colormap gray; axis image; title(method_names{meth});
        end
    end

    x = x0;
    total_iter = 0;
    tic;
    for k = 1:K
        disp("Running method: SOR for size: " + m + " run: " + k)
        A = FormMatrix(x, alpha);
        b = u0;
        [x, iter] = SOR(omega_opt(s), A, b, x, maxiter, tol);
        total_iter = total_iter + iter;
    end
    elapsed = toc;
    results_time(s, end) = elapsed;
    results_iter(s, end) = total_iter;
    if s == 1
        figure; imagesc(reshape(x, m, m)); colormap gray; axis image; title('SOR');
    end
end

% results
disp('Execution times (seconds):');
disp(array2table(results_time, 'VariableNames', [method_names, {'SOR'}], 'RowNames', cellstr(num2str(sizes'))));
disp('Total iterations:');
disp(array2table(results_iter, 'VariableNames', [method_names, {'SOR'}], 'RowNames', cellstr(num2str(sizes'))));

profile viewer;
