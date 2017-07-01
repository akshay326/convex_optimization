%% simple_portfolio_data
rand('state', 5);
randn('state', 5);
n=20;
pbar = ones(n,1)*.03+[rand(n-1,1); 0]*.12;
S = randn(n,n);
S = S'*S;
S = S/max(abs(diag(S)))*.2;
S(:,n) = zeros(n,1);
S(n,:) = zeros(n,1)';
x = ones(n,1)/n;

% optimizing
cvx_begin
    variable x(n)
    minimize(x'*S*x)
    subject to
        pbar'*x >= 0
        1 >= ones(1,n)*x >= 1
        x >= 0
cvx_end