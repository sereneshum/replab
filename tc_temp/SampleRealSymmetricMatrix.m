function T = SampleRealSymmetricMatrix(n)
    %
    % Generates a real symmetric matrix with measure invariant 
    % under orthogonal transformations.
    %
    % Reference:
    %     http://staff.math.su.se/shapiro/UIUC/random_matrices.pdf
    %     Section 4. Gaussian Ensembles, second definition page 20
    %

    % initialization of n x n matrix
    T = zeros(n, n);

    for i = 1:n
        % diagonal elements set to iid normal variables
        % with variance scaled up by a factor of sqrt(2)
        T(i, i) = randn * sqrt(2);

        % off-diagonal elements iid standard normal variables
        for j = 1:i-1 % if 1 = 0, then this loop does nothing
            T(i, j) = randn;
            T(j, i) = T(i, j);
        end
    end
end
