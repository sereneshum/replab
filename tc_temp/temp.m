disp('Hello world')
x=sqrt(16)
a=SampleGinibreEnsemble(3)
b=SampleUnitaryMatrix(10)
% ctranspose(b)*b
% eigvalues=eig(b)
% plot(eigvalues,'*')
% for i = 1:0
%     disp(i);
% end
T = SampleRealSymmetricMatrix(5)
v = SampleComplexVector(5)
[a,b,q]=UnstableLanczosIteration(T, v, 10)

A = SampleHermitianMatrix(10)
[a,b,V]=LanczosTridiagonalization(kron(A, eye(5)))

A = eye(10)
[a,b,V]=LanczosTridiagonalization(A)
