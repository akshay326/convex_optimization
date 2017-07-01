% A = randn(2,2);
% [U,S,V] = svd2x2(A);

cvx_begin
    variable x(2)
    minimize quad_form(x,V)
    subject to
        1 <= x <=2;
cvx_end