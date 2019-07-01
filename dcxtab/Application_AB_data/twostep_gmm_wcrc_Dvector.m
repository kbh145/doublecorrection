function [b2,sr,sw,s0,b2c,src,swc,s0c,b1,s1r,s10,J,pv] = twostep_gmm_wcrc_Dvector(y,x,z,mem)

% Calculates one-step and two-step GMM estimator for linear model
% y = X*b + e with instrument matrix Z
% Calculates covariance matrix and standard errors robust to
% misspecification and clusters
% Also calculates Windmeijer corrected standard errors
%
% Under random sampling, one-step GMM is 2SLS
% Under cluster sampling, one-step GMM is Arellano-Bond one-step GMM (ZHZ)
% Inputs:
%	y       nx1 vector of observations
%	x       nxk matrix of regressors
%	z       nxl matrix of instruments, l>=k (includes exogenous components of x)
%   mem     nx1 vector of cluster membership 
%           set mem=(1:n)' for iid
%           set mem=kron((1:N)',ones(Tbar,1)) for balanced panel where
%           N is the number of cross-section units and Tbar is the time
%           series length used for estimation

% Outputs:
%	b2      kx1 vector of coefficient estimates with uncentered w
%	sr      kx1 vector of misspecification-and-cluster-and-heteroskedasticity
%           robust asymptotic standard errors
%	sw      kx1 vector of Windmeijer-corrected cluster-and-heteroskedasticity
%           robust asymptotic standard errors
%	s0      kx1 vector of classic cluster-and-heteroskedasticity robust
%           asymptotic standard errors
%
%	b2c     kx1 vector of coefficient estimates with centered weight matrix
%	src      kx1 vector of misspecification-and-cluster-and-heteroskedasticity
%           robust asymptotic standard errors with centered weight matrix
%	swc     kx1 vector of Windmeijer-corrected cluster-and-heteroskedasticity
%           robust asymptotic standard errors with centered weight matrix
%	s0c     kx1 vector of classic cluster-and-heteroskedasticity robust
%           asymptotic standard errors with centered weight matrix
%
%	J       J-statistic for overidentifying restrictions 
%	pv      Asymptotic chi-square p-value for J-statistic


n = size(y,1);
k = size(x,2);
l = size(z,2);
G = length(unique(mem));

idx = repmat(mem,1,G)==kron(ones(n,1),unique(mem)');
ng = sum(idx)';

zx = z'*x;
zy = z'*y;

if n == G
    w0 = (z'*z)/n;
else
    w0 = zeros(l,l);

    for g = 1:G
        zg = z(idx(:,g),:);
    
        temp_H = [zeros(1,ng(g)-1) ; diag(ones(ng(g)-1,1)) ];
        temp_H = [temp_H zeros(ng(g),1)]* -1;
        H = temp_H + temp_H' + diag(ones(ng(g),1))*2;
    
        w0 = w0 + zg'*H*zg;
    end

    w0 = w0/n;
end

b1 = (zx'/w0*zx)\(zx'/w0*zy);

cn = sum(ng.^2)/n;

e1 = y - x*b1;
ze1 = z.*(e1*ones(1,l));
mu1 = mean(ze1)';

w = zeros(l,l);

if n == G
   w = (ze1'*ze1)/n;
   wc = w - cn*(mu1*mu1');
else
   for g = 1:G
       zg = z(idx(:,g),:);
       eg = e1(idx(:,g));
       zeg = zg'*eg;
       
       w = w + zeg*zeg';
   end
   w = w/n;
end

wc = w - cn*(mu1*mu1');

b2 = (zx'/w*zx)\(zx'/w*zy);
b2c = (zx'/wc*zx)\(zx'/wc*zy);

e2 = y - x*b2;
e2c = y - x*b2c;
ze2 = z.*(e2*ones(1,l));
mu2 = mean(ze2)';
ze2c = z.*(e2c*ones(1,l));
mu2c = mean(ze2c)';

if l>k
  J = (mu2'/w*mu2)*n;
  pv = chi2cdf(J,l-k,'upper');
%   Jc = (mu'/wc*mu)*n;
%   pvc = chi2cdf(Jc,l-k,'upper');
else
  J = 0;
  pv = 1;
%   Jc = 0;
%   pvc = 1;
end

H0 = (1/n^2)*zx'/w0*zx;
H1 = (1/n^2)*zx'/w*zx;
H1c = (1/n^2)*zx'/wc*zx;

if n == G
    D = (zx'/w*zx)\zx'/w*(2/n)*(z'*(x.*repmat(ezwze,1,k)));
    Dc = (zx'/wc*zx)\zx'/wc*(2/n)*(z'*(x.*repmat(ezwcze,1,k))-(1/n^2)*zx'*mu1/w*mu2c-(1/n^2)*mu1*mu2c'/wc*zx);
    Psi1 = (1/n)*(z.*repmat(e1,1,l))/w0*(z'*x)+(1/n)*repmat(((e1'*z)/w0*z')',1,k).*x-(1/n^2)*repmat(((e1'*z)/w0*z')',1,k).*(z/w0*z'*x);
    Psi2 = (1/n)*(z.*repmat(e2,1,l))/w*(z'*x)+(1/n)*repmat(((e2'*z)/w*z')',1,k).*x-(1/n^2)*repmat(((e2'*z)/w*z')',1,k).*((z.*repmat(e1.^2,1,l))/w*z'*x);
    Psi2c = (1/n)*(z.*repmat(e2c,1,l))/wc*(z'*x)+(1/n)*repmat(((e2c'*z)/wc*z')',1,k).*x-(1/n^2)*repmat(((e2c'*z)/wc*z')',1,k).*((z.*repmat(e1.^2,1,l))/wc*z'*x);
else
    D_part = zeros(k,k);
    Dc_part = zeros(k,k);
    Psi1 = zeros(G,k);
    Psi2 = zeros(G,k);
    Psi2c = zeros(G,k);
    for g = 1:G
        zg = z(idx(:,g),:);
        e1g = e1(idx(:,g));
        e2g = e2(idx(:,g));
        e2cg = e2c(idx(:,g));
        xg = x(idx(:,g),:);
        
        zge1 = zg.*(e1g*ones(1,l));
        mug1 = mean(zge1)';
        ng = length(zg(:,1));
        muzxg = zg'*xg/ng;
        
        D_part = D_part + kron((e2'*z)/w,zx'/w)* (kron(zg'*e1g,zg'*xg) + kron(zg'*xg,zg'*e1g));
        Dc_part = Dc_part + kron((e2c'*z)/wc,zx'/wc)* (kron(zg'*e1g,zg'*xg) + kron(zg'*xg,zg'*e1g));
        %D_part2 = D_part + (zg'*e1g)*(e2'*z)/w*(zg'*xg) + (zg'*xg)*((e2'*z)/w*(zg'*e1g))
        %Dc_part2 = Dc_part + (zg'*e1g)*(e2c'*z)/wc*(zg'*xg) + (zg'*xg)*((e2c'*z)/wc*(zg'*e1g))

        
        temp_H = [zeros(1,ng-1) ; diag(ones(ng-1,1)) ]; temp_H = [temp_H zeros(ng,1)]* -1;
        H = temp_H + temp_H' + diag(ones(ng,1))*2;

        Psi1(g,:) = ((1/n)*(x'*z)/w0*(zg'*e1g)+(1/n)*(xg'*zg)/w0*(z'*e1)-(1/n^2)*(x'*z)/w0*(zg'*H*zg)/w0*(z'*e1))';
        Psi2(g,:) = ((1/n)*(x'*z)/w*(zg'*e2g)+(1/n)*(xg'*zg)/w*(z'*e2)-(1/n^2)*(x'*z)/w*(zg'*e1g)*(e1g'*zg)/w*(z'*e2))';
        Psi2c(g,:) = ((1/n)*(x'*z)/wc*(zg'*e2cg)+(1/n)*(xg'*zg)/wc*(z'*e2c)-(1/n^2)*(x'*z)/wc*(((zg'*e1g)-mu1)*((e1g'*zg)-mu1'))/wc*(z'*e2c))';
    end
    
    D_part = D_part/n;
    D = n^2*H1\D_part;
    %D = n^2*H1\zx'/w*D_part;
    Dc_part = Dc_part/n;
    Dc_part = Dc_part - kron((e2c'*z)/wc,zx'/wc)*(cn*kron(mu1,zx/n)-cn*kron(zx/n,mu1));
%    Dc_part = Dc_part - cn*(zx/n)*mu1'/wc*(z'*e2c) - cn*mu1*(e2c'*z)/wc*(zx/n);
    Dc = n^2*H1c\Dc_part;
    %Dc = n^2*H1c\zx'/wc*Dc_part;
end

Om1 = (Psi1'*Psi1)/n;
Om2 = (Psi2'*Psi2)/n;
Om2c = (Psi2c'*Psi2c)/n;
Om12 = (Psi1'*Psi2)/n;
Om12c = (Psi1'*Psi2c)/n;

Vtil2 = H1\Om2/H1';
Vtil1 = H0\Om1/H0';
Ctil = H0\Om12/H1';

Vtil2c = H1c\Om2c/H1c';
Ctilc = H0\Om12c/H1c';

s1r = sqrt(diag(Vtil1/n));

Vr = Vtil2 + D*Ctil + Ctil'*D' + D*Vtil1*D';
sr = sqrt(diag(Vr/n));

Vrc = Vtil2c + Dc*Ctilc + Ctilc'*Dc' + Dc*Vtil1*Dc';
src = sqrt(diag(Vrc/n));

Q = -zx/n;
V0 = inv(Q'/w*Q);
s0 = sqrt(diag(V0/n));
V01 = H0\(Q'/w0*w/w0*Q)/H0';
V01c = H0\(Q'/w0*wc/w0*Q)/H0';

s10 = sqrt(diag(V01/n));

V0c = inv(Q'/wc*Q);
s0c = sqrt(diag(V0c/n));

Vw = V0 + D*V0 + V0*D' + D*V01*D';
sw = sqrt(diag(Vw/n));

Vwc = V0c + Dc*V0c + V0c*Dc' + Dc*V01c*Dc';
swc = sqrt(diag(Vwc/n));

