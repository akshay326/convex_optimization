%bnds = randn(n,2);
%l = min( bnds, [], 2 );
%u = max( bnds, [], 2 );

%cvx_begin
%    variable x(n)
%    minimize(norm(A*x - b))
%    subject to
%    l <= x <= u
%cvx_end

cvx_begin
    variable x(1)
    variable y(1)
    minimize(x*x + 9*y*y)
    subject to
        0 <= x
        0 <= y
        2*x + y >= 1
        x + 3*y >= 1
cvx_end