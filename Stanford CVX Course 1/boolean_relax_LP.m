rng(0,'v5uniform');
n=100;
m=300;
A=rand(m,n);
b=A*ones(n,1)/2;
c=-rand(n,1);

% finding an Lower Bound
cvx_begin
    variable x(n)
    minimize(c'*x)
    subject to
        A*x <= b
        0 <= x <= 1
cvx_end

% guessing via plotting
figure % opens new figure window
L = ones(1,100);
C = ones(1,100);
for i = 1:100
   y = x > (i/100);
   L(1,i) = c'*y - c'*x;
   C(1,i) = c'*y
end
plot(0:0.01:0.99,L,'--',0:0.01:0.99,C)