output_precision(3)
disp('Lanczos Testing Branch!')
a=SampleGinibreEnsemble(3)
b=SampleUnitaryMatrix(10)
b'*b % check if b is unitary 
eigvalues=eig(b)
% plot(eigvalues,'*')
% T = SampleRealSymmetricMatrix(5)
% v = SampleComplexVector(5)
% [a,b,q]=UnstableLanczosIteration(T, v, 10)

A = kron(SampleHermitianMatrix(7), eye(20));
[a,b,V]=LanczosTridiagonalization(A);

m=length(a);
T = zeros(m,m);
for i = 1:m
    T(i,i) = a(i);
end
for i = 1:m-1
    T(i,i+1) = b(i);
    T(i+1,i) = b(i);
end
disp('dimensions of hermitian A matrix:')
size(A)
disp('dimensions of tridiagonal T matrix:')
size(T)
disp('norm error ||V''AV - T||')
norm(V'*A*V - T)
disp('eigenvalues of T, and A')
(round(eig(T).*1e3)./1e3).'
unique(round(eig(A).*1e3)./1e3).'

% A = eye(10)
% [a,b,V]=LanczosTridiagonalization(A)
