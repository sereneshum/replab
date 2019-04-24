function [ a, b, V ] = LanczosTridiagonalization( A, esp )

    %   Computes the Lanczos Tridiagonalization of A using complete
    %   reorthonormalization.
    %
    %   The value m is the number of distinct eigenvalues of A and is
    %   determined computationally.
    %
    %   Bulk equation (1 < j < m):
    %       A v(j) = b(j) v(j+1) + a(j) v(j) + b(j-1) v(j-1) where v(j) = V(:,j)
    %
    %   Returns:
    %       a :: [1 x m] row-vector of diagonal values of T
    %       b :: [1 x m] row-vector of off-diagonal values of T (always starts with zero)
    %       V :: [n x m] isometry matrix which satisfies V' A V = T

    % Parsing input parameters
    if nargin < 1
        error('LanczosTridiagonalization requires input matrix A that is Hermitian');
    end
    if nargin < 2
        esp = 1e-10;
    end

    % since m <= n always, allocate enough space for the worst case scenario
    % truncation of a,b,Q will occur at the end
    m = -1;
    n = length(A);
    a = zeros(1,n);
    b = zeros(1,n);
    V = zeros(n,n);

    % abbreviated initial iteration step
    V(:,1) = ones(n,1)/sqrt(n); % following C. Guo and S. Qiao, July 2003, Section 5

    % main Lanczos iteration
    for j = 1:n-1
        % computes V(:j+1), a(j) and b(j)
        % using all previous values

        % v is a [nx1] column vector
        % which is supposed to satisfy
        % 1) v is orthogonal to all V(:,k) for all k < j+1
        % 2) v is in the Krylov subspace of A w.r.t. V(:,1)
        v = A * V(:,j);
        a(j) = real(V(:,j)' * v); % since A is Hermitian, a is always real
        v = v - a(j) * V(:,j);

        % TODO: this step seems pointless because of the reorthogonalization
        if j > 1
            v = v - b(j-1) * V(:,j-1);
        end

        % in order to ensure orthonormality of columns of V,
        % reorthogonalize with respect to previous vectors
        for i = 1:j-1
            v = v - ( V(:,j-i)' * v ) * V(:,j-i);
        end
        % and then normalize if its a substantial vector
        b(j) = norm(v);
        if b(j) < esp
            % if the new vector is smaller than some tolerance,
            % then the Krylov subspace has been spanned, stop
            m = j;
            break
        end
        V(:,j+1) = v/b(j);
    end

    % if m < 0, then it was never changed from m = -1
    % therefore, m = n and no truncation of a, b, or V is needed
    if m > 0
        % truncate a, b, V according to the value of m
        a = a(1:m);
        b = b(1:m-1); % b(m) is defined at this point, but close to zero and not needed
        V = V(:,1:m);
    end
end
