function D2 = distmetric_binomial_probability( X1, X2 )
% finds distance of 2 vectors using binomial probability
% The returned values is basically a p value.
% X1 is an 1 x n row vector
% X2 is a m2 x n matrix
% D2 is a m2 x 1 vector representing distances

% Turn feature into binary
coverage_thresh = 11;
A = X1>coverage_thresh;
X2 = X2>coverage_thresh;

numrows = size(X2,1);
D2 = zeros(numrows,1);
tot = length(A);

for row = 1:numrows
    B = X2(row,:);
    % create 2x2 count matrix
    tt = sum(A & B);
    tn = sum(A & ~B);
    nt = sum(~A & B);
    nn = sum(~(A | B));
    t1 = tt+tn;
    n1 = nt+nn;
    t2 = tt+nt;
    n2 = tn+nn;
    assert(tot == t1+n1);
    assert(tot == t2+n2);
    D2(row) = factorial(n1)/factorial(tot)*factorial(n2)/factorial(nn)/...
        factorial(nt)*factorial(t1)/factorial(tn)*factorial(t2)/factorial(tt);
end

end

