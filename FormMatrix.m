function A = FormMatrix(u, alpha)
    beta = 1e-6;
    n = length(u);
    m = sqrt(n);
    h = 1 / (m + 1);

    U = reshape(u, m, m);

    maxEntries = 5 * n;
    I = zeros(maxEntries,1);
    J = zeros(maxEntries,1);
    V = zeros(maxEntries,1);
    idx_count = 0;

    for j = 1:m
        for i = 1:m
            idx = (j-1)*m + i;

            u_ij = U(i,j);

            % Precompute neighbors
            if i > 1
                u_im1j = U(i-1,j);
            else
                u_im1j = 0;
            end
            if i < m
                u_ip1j = U(i+1,j);
            else
                u_ip1j = 0;
            end
            if j > 1
                u_ijm1 = U(i,j-1);
            else
                u_ijm1 = 0;
            end
            if j < m
                u_ijp1 = U(i,j+1);
            else
                u_ijp1 = 0;
            end
            if i < m && j > 1
                u_ip1jm1 = U(i+1,j-1);
            else
                u_ip1jm1 = 0;
            end
            if i > 1 && j < m
                u_im1jp1 = U(i-1,j+1);
            else
                u_im1jp1 = 0;
            end

            % differences and squares to prevent recomputation
            dux_w = (u_ij - u_im1j)/h;
            dux_w2 = dux_w^2;
            dux_e = (u_ip1j - u_ij)/h;
            dux_e2 = dux_e^2;
            duy_s = (u_ij - u_ijm1)/h;
            duy_s2 = duy_s^2;
            duy_em1 = (u_ip1j - u_ip1jm1)/h;
            duy_em1_2 = duy_em1^2;
            duy_ep1 = (u_ijp1 - u_ij)/h;
            duy_ep1_2 = duy_ep1^2;
            duy_n = (u_ijp1 - u_im1jp1)/h;
            duy_n2 = duy_n^2;
            du_ip1jm1 = (u_ip1jm1 - u_ijm1)/h;
            du_ip1jm1_2 = du_ip1jm1^2;

            % AW
            if i > 1
                t1 = sqrt(dux_w2 + duy_s2 + beta);
                if j < m
                    t2 = sqrt(dux_w2 + ((u_ijp1 - u_im1j)/h)^2 + beta);
                else
                    t2 = sqrt(dux_w2 + beta);
                end
                AW = -alpha/(h^2) * 0.5 * (1/t1 + 1/t2);
            else
                AW = 0;
            end

            % AE
            if i < m
                t1 = sqrt(dux_e2 + duy_em1_2 + beta);
                if j < m
                    t2 = sqrt(dux_e2 + duy_ep1_2 + beta);
                else
                    t2 = sqrt(dux_e2 + beta);
                end
                AE = -alpha/(h^2) * 0.5 * (1/t1 + 1/t2);
            else
                AE = 0;
            end

            % AS
            if j > 1
                t1 = sqrt(dux_w2 + duy_s2 + beta);
                if i < m
                    t2 = sqrt(du_ip1jm1_2 + duy_s2 + beta);
                else
                    t2 = sqrt(duy_s2 + beta);
                end
                AS = -alpha/(h^2) * 0.5 * (1/t1 + 1/t2);
            else
                AS = 0;
            end

            % AN
            if j < m
                t1 = sqrt(dux_e2 + duy_ep1_2 + beta);
                if i > 1
                    t2 = sqrt(duy_n2 + duy_ep1_2 + beta);
                else
                    t2 = sqrt(duy_ep1_2 + beta);
                end
                AN = -alpha/(h^2) * 0.5 * (1/t1 + 1/t2);
            else
                AN = 0;
            end

            AC = -(AW + AE + AS + AN) + 1;

            % Diagonal
            idx_count = idx_count + 1;
            I(idx_count) = idx;
            J(idx_count) = idx;
            V(idx_count) = AC;

            % West 
            if i > 1
                idx_count = idx_count + 1;
                I(idx_count) = idx;
                J(idx_count) = idx - 1;
                V(idx_count) = AW;
            end

            % East 
            if i < m
                idx_count = idx_count + 1;
                I(idx_count) = idx;
                J(idx_count) = idx + 1;
                V(idx_count) = AE;
            end

            % South 
            if j > 1
                idx_count = idx_count + 1;
                I(idx_count) = idx;
                J(idx_count) = idx - m;
                V(idx_count) = AS;
            end

            % North 
            if j < m
                idx_count = idx_count + 1;
                I(idx_count) = idx;
                J(idx_count) = idx + m;
                V(idx_count) = AN;
            end
        end
    end

    I = I(1:idx_count);
    J = J(1:idx_count);
    V = V(1:idx_count);

    A = sparse(I, J, V, n, n);
end