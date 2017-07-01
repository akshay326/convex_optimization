% Target Expresion : x*x + 1
cvx_begin
    variable x(1);
    minimize x*x + 1
    subject to
        (x-2)*(x-4) <= 0;
cvx_end

% Plotting Part
% x = -1:0.1:5;
% f0 = x.^2 + 1;
% f1 = (x-2).*(x-4);
% L1 = x.^2 + 1 + (x-2).*(x-4);
% L2 = x.^2 + 1 + 1.5*(x-2).*(x-4);
% L3 = x.^2 + 1 + 0.5*(x-2).*(x-4);
% plot(x,f0,x,f1,x,L1,x,L2,x,L3);
% grid on;