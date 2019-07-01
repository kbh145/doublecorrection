function [b,s,sc,sw,s0,J,pv,iter] = iterated_gmm_wcrc(y,x,z,mem)

% Calculates iterated GMM estimator for linear model   y = X*b + e with instrument matrix Z
% Calculates covariance matrix and standard errors robust to
% misspecification and clusters
% Also calculates Windmeijer corrected standard errors

% Inputs:
%	y       nx1 vector of observations
%	x       nxk matrix of regressors
%	z       nxl matrix of instruments, l>=k (includes exogenous components of x)
%   mem     nx1 vector of cluster membership set mem=(1:n)' for iid

% Outputs:
%	b       kx1 vector of coefficient estimates
%	s       kx1 vector of misspecification-and-cluster-and-heteroskedasticity
%           robust asymptotic standard errors
%	sc      kx1 vector of misspecification-and-cluster-and-heteroskedasticity
%           robust asymptotic standard errors with centered weight matrix
%	sw      kx1 vector of Windmeijer-corrected cluster-and-heteroskedasticity
%           robust asymptotic standard errors
%	s0      kx1 vector of classic cluster-and-heteroskedasticity robust
%           asymptotic standard errors
%	J       J-statistic for overidentifying restrictions 
%	pv      Asymptotic chi-square p-value for J-statistic
%
% Output variables = NaN if the iteration reaches maxit

tolerance = 1e-5;
maxit = 1e+3;

n = size(y,1);
k = size(x,2);
l = size(z,2);
G = length(unique(mem));
zx = z'*x;
zy = z'*y;
w = z'*z;
b1 = (zx'/w*zx)\(zx'/w*zy);

idx = repmat(mem,1,G)==kron(ones(n,1),unique(mem)');

ng = sum(idx)';
cn = sum(ng.^2)/n;

for iter = 1:maxit
   e = y - x*b1;
   w = zeros(l,l);
   
   if n == G
       ze = z.*repmat(e,1,l);
       w = (ze'*ze)/n;
   else
       for g = 1:G
           zg = z(idx(:,g),:);
           eg = e(idx(:,g));
           zeg = zg'*eg;
           w = w + zeg*zeg';
       end
       w = w/n;
   end
   
   b = (zx'/w*zx)\(zx'/w*zy);
   db = b - b1;
   if norm(db) < tolerance
       break
   end   
   b1 = b;
   
   if iter == maxit
       b = NaN;
       s = NaN;
       V = NaN;
       sc = NaN;
       Vc = NaN;
       sw = NaN;
       Vw = NaN;
       s0 = NaN;
       V0 = NaN;
       J = NaN;
       pv = NaN;
       Jc = NaN;
       pvc = NaN;
       return
   end
end

e = y - x*b;
ze = z.*(e*ones(1,l));
mu = mean(ze)';

wc = w - cn*(mu*mu');

if l>k
  J = (mu'/w*mu)*n;
  pv = chi2cdf(J,l-k,'upper');
  Jc = (mu'/wc*mu)*n;
  pvc = chi2cdf(Jc,l-k,'upper');
else
  J = 0;
  Jc = 0;
  pv = 1;
  pvc = 1;
end

if n == G
    ezwze = e.*(z/w*z'*e);
    H = (1/n^2)*(x'*z)/w*(z'*x)-(2/n^3)*(x'*z)/w*(z'*(x.*repmat(ezwze,1,k)));
    Hc = (1+cn*Jc/n)*H;
    Psi = -(1/n)*(z.*repmat(e,1,l))/w*(z'*x)-(1/n)*repmat(((e'*z)/w*z')',1,k).*x+(1/n^2)*repmat(((e'*z)/w*z')',1,k).*((z.*repmat(e.^2,1,l))/w*z'*x);
    Psic = (1+cn*Jc/n)*Psi + (1/n)*(cn*Jc/n)*(z.*repmat(e,1,l))/w*(z'*x);
else
    Hpart = zeros(l,k);
    Psi = zeros(G,k);
    Psic = zeros(G,k);
    for g = 1:G
        zg = z(idx(:,g),:);
        eg = e(idx(:,g));
        xg = x(idx(:,g),:);
        Hpart = Hpart + (zg'*eg)*(e'*z)/w*(zg'*xg) + (zg'*xg)*((e'*z)/w*(zg'*eg));

        Psi(g,:) = (-(1/n)*(x'*z)/w*(zg'*eg)-(1/n)*(xg'*zg)/w*(z'*e)+(1/n^2)*(x'*z)/w*(zg'*eg)*(eg'*zg)/w*(z'*e))';
        Psic(g,:) = (1+cn*Jc/n)*Psi(g,:) + ((1/n)*(cn*Jc/n)*(x'*z)/w*(zg'*eg))';
    end
    H = (1/n^2)*(x'*z)/w*(z'*x)-(1/n^3)*(x'*z)/w*Hpart;
    Hc = (1+cn*Jc/n)*H;

end

Om = (Psi'*Psi)/n;
Omc = (Psic'*Psic)/n;

V = H\Om/H';
s = sqrt(diag(V/n));

Vc = Hc\Omc/Hc';
sc = sqrt(diag(Vc/n));

Q = -zx/n;
V0 = inv(Q'/w*Q);
s0 = sqrt(diag(V0/n));

Vw = H\(Q'/w*Q)/H';
sw = sqrt(diag(Vw/n));

