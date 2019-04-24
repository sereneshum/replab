function [ a, b, q ] = UnstableLanczosIteration( T, f, m )

    a = zeros(1,1+m);
    b = zeros(1,m);
    q = zeros(length(f),m+1);

    q(:,2) = f/norm(f);
    for i = 2:m
        v = T*q(:,i);
        a(i) = real(ctranspose(q(:,i))*v); % a should always be real
        v = v - b(i-1)*q(:,i-1) - a(i)*q(:,i);
        b(i) = norm(v);
        if b(i) ~= 0
            q(:,i+1) = v/b(i);
        else
            % q(:,i+1) is already initialized to zero
        end
    end

    % The initial values of a, b, and q are all zero
    a = a(2:end);
    b = b(2:end);
    q = q(:,2:end);
end
