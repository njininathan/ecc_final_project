function H = obtain_H(base_matrix)
    
    [M, N] = size(base_matrix);
    Q = 96;
    H = zeros(Q * M, Q * N);
    I = eye(Q);
    
    for i = 1:M
        for j = 1:N
            row_range = (i - 1) * Q + 1 : i * Q;
            col_range = (j - 1) * Q + 1 : j * Q;
            if base_matrix(i, j) >=0
                H(row_range, col_range) = circshift(I, base_matrix(i, j), 2); % Only apply shift if value is non-negative
            end
        end
    end

end

