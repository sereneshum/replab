function T = SampleHermitianMatrix(n)
    %
    % Generates a Hermitian matrix with measure invariant 
    % under unitary transformations.
    %
    % Reference:
    %     http://staff.math.su.se/shapiro/UIUC/random_matrices.pdf
    %     Section 4. Gaussian Ensembles, first definition page 20
    %

    % initialization of n x n matrix
    T = zeros(n, n);

    for i = 1:n
        % diagonal elements set to iid standard normal variables
        T(i, i) = randn;

        % off-diagonal elements iid standard normal variables
        for j = 1:i-1 % if 1 = 0, then this loop does nothing
            T(i, j) = complex(randn, randn)/sqrt(2);
            T(j, i) = conj(T(i, j));
        end
    end
end
