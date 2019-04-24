function Q = SampleUnitaryMatrix(n)
    [Q,R] = qr(SampleGinibreEnsemble(n));
end
