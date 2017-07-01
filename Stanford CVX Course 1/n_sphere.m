x = rand(2,3)
cvx_begin
    maximize(norms(norms(x).^-1))
    subject to
        1 <= x*x' <= 1;
        x(2)*x(1)' <= 0.5;
        x(3)*x(1)' <= 0.5;
        x(2)*x(3)' <= 0.5;
cvx_end
x
% x = random_Nsphere(3);
