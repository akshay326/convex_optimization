% Minimax/Quadratic Fractions using Bisection Method
k=201;
t=(-3:6/(k-1):3)';
y=exp(t);
Tpowers=[ones(k,1) t t.^2];

u=exp(3);
l=0; % initial upper and lower bounds
bisection_tol=1e-3; % bisection tolerance

while u-l>= bisection_tol
    gamma=(l+u)/2;
    
    cvx_begin % solve the feasibility problem
    cvx_quiet(true);
    variable a(3);
    variable b(2);
    subject to
    abs(Tpowers*a-y.*(Tpowers*[1;b])) <= gamma*Tpowers*[1;b];
    cvx_end
    
    if strcmp(cvx_status,'Solved')
        u=gamma;
        a_opt=a;
        b_opt=b;
        objval_opt=gamma;
    else
        l=gamma;    
    end
end

% Mapping it
y_fit=Tpowers*a_opt./(Tpowers*[1;b_opt]);
figure(1);
plot(t,y,'b', t,y_fit,'r+');
xlabel('t');
ylabel('y');
figure(2);
plot(t, y_fit-y);
xlabel('t');
ylabel('err');