% Target Expresion : x*x + 2*y*y - x*y - x
P = [1 , -1, -1;0,2,0;0,0,1];

% For Pertubation Analysis, uncomment issue pert
% Issue pert
% d = [-0.1,0,0.1];

cvx_begin
    variable x(3);
    dual variables lambda{3};
    minimize quad_form(X,P)
    subject to
        lambda{1} : x(1) + 2*x(2) <= -2+d(i);
        lambda{2} : x(1) - 4*x(2) <= -3+d(j);
        lambda{3} : x(1) + x(2) >= -5;
        x(3) == 1;
cvx_end

% Issue Pert
% for i = 1:3
%     for j = 1:3
% cvx_begin
%     variable x(3);
%     dual variables lambda{3};
%     minimize quad_form(X,P)
%     subject to
%         lambda{1} : x(1) + 2*x(2) <= -2+d(i);
%         lambda{2} : x(1) - 4*x(2) <= -3+d(j);
%         lambda{3} : x(1) + x(2) >= -5;
%         x(3) == 1;
% cvx_end
%     input('contniue') % request input
%     end
% end