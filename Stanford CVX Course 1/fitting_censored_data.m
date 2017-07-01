% data for censored fitting problem.
randn('state',0);

n = 20;  % dimension of x's
M = 25;  % number of non-censored data points
K = 100; % total number of points
c_true = randn(n,1);
X = randn(n,K);
y = X'*c_true + 0.1*(sqrt(n))*randn(K,1);

% Reorder measurements, then censor
[y, sort_ind] = sort(y);
X = X(:,sort_ind);

% Exclude Censored Data 
y = y(1:M);
X = X(:,1:M);

% Include Censored Data 
% D = (y(M)+y(M+1))/2;
% y(M+1:K) = D;

cvx_begin
    variable c(n)
    minimize norm(y - X'*c)
cvx_end