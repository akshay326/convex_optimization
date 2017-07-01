mu1 = 8; mu2 = 20;
sigma1 = 6; sigma2 = 17.5;
rho = -0.25;

% assuming jointly gaussian distribution
mu = mu1 + mu2;
sigma = sqrt(sigma1^2 + sigma2^2 + 2*rho*sigma1*sigma2);
ploss = normcdf(0, mu, sigma); % gaussian probability of loss
n = 100;
rmin = -30; rmax = 70;

% discretize outcomes of R1 and R2
r = linspace(rmin,rmax,n)';

% marginal distributions
p1 = exp(-(r-mu1).^2/(2*sigma1^2)); p1 = p1/sum(p1);
p2 = exp(-(r-mu2).^2/(2*sigma2^2)); p2 = p2/sum(p2);

% form mask of region where R1 + R2 <= 0
r1p = r*ones(1,n); r2p = ones(n,1)*r';
loss_mask = (r1p + r2p <= 0)';


cvx_begin
    variable P(n,n)
    maximize (sum(sum(P(loss_mask))))
    subject to
        P >= 0;
        sum(P ,2) == p1;
        sum(P',2) == p2;
        (r - mu1)'*P*(r - mu2) == rho*sigma1*sigma2;
cvx_end


pmax = cvx_optval; % worst case probability of loss
pmax
ploss
P = full(P);

% Mesh
figure(1); mesh(r1p, r2p, P');
xlabel('R1'); ylabel('R2'); zlabel('density');

% Contours
xlim([rmin rmax]); ylim([rmin rmax]);
print -depsc 'loss_bounds_mesh.eps'
figure(2); contour(r1p, r2p, P');
xlabel('R1'); ylabel('R2'); grid on;
xlim([rmin rmax]); ylim([rmin rmax]);
print -depsc 'loss_bounds_cont.eps'