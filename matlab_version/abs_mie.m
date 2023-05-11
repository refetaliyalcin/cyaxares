function [Qext,w0,g,s1,s2,P11,QI,S_FACTOR,SCALE_FACTOR] = abs_mie(radius,lambda,n_medium,k_medium,n_pigment,k_pigment,nang)

m1 = n_medium + 1i*k_medium;
m2 = n_pigment + 1i*k_pigment;
% Area=pi*radius^2;
x0=2*pi*radius/lambda; %size parameter at free-space

mr = m2/m1;
x = m1*x0;
y = m2*x0;
nmax = ceil(abs(x) + 4.05*abs(x)^0.3333 + 8.0);
n=(1:nmax);
n_wide=0:nmax+1;

nu = (n_wide+0.5)';

besselj_x=besselj(nu,x);
besselj_x_minus_1=besselj_x;
besselj_x_minus_1(end-1:end)=[];
besselj_x_plus_1=besselj_x;
besselj_x_plus_1(1:2)=[];
besselj_x_middle=besselj_x;
besselj_x_middle(end)=[];
besselj_x_middle(1)=[];

besselj_y=besselj(nu,y);
besselj_y_minus_1=besselj_y;
besselj_y_minus_1(end-1:end)=[];
besselj_y_plus_1=besselj_y;
besselj_y_plus_1(1:2)=[];
besselj_y_middle=besselj_y;
besselj_y_middle(end)=[];
besselj_y_middle(1)=[];

bessely_x=bessely(nu,x);
bessely_x_minus_1=bessely_x;
bessely_x_minus_1(end-1:end)=[];
bessely_x_plus_1=bessely_x;
bessely_x_plus_1(1:2)=[];
bessely_x_middle=bessely_x;
bessely_x_middle(end)=[];
bessely_x_middle(1)=[];


psi=sqrt(0.5*pi*x)*besselj_x_middle;
dpsi_chain_1=sqrt(0.5*pi*x)*0.5*(besselj_x_minus_1 - besselj_x_plus_1);
dpsi_chain_2=0.25*pi*besselj_x_middle/sqrt(0.5*pi*x);
psid=dpsi_chain_1 + dpsi_chain_2;

psi2=sqrt(0.5*pi*y)*besselj_y_middle;
dpsi_chain_1=sqrt(0.5*pi*y)*0.5*(besselj_y_minus_1 - besselj_y_plus_1);
dpsi_chain_2=0.25*pi*besselj_y_middle/sqrt(0.5*pi*y);
psi2d=dpsi_chain_1 + dpsi_chain_2;


xi=sqrt(0.5*pi*x)*bessely_x_middle;
dxi_chain_1=sqrt(0.5*pi*x)*0.5*(bessely_x_minus_1 - bessely_x_plus_1);
dxi_chain_2=0.25*pi*bessely_x_middle/sqrt(0.5*pi*x);
xid=dxi_chain_1 + dxi_chain_2;


xi = psi +  xi*1j;
xid = psid + xid*1j;

ct = psi.*psi2d;
ct2 = psi2.*psid;
% ct3 = mr*(psi.*xid - xi.*psid);
ct4 = (mr*psi2.*xid - xi.*psi2d);
ct5 = psi2.*xid - mr*xi.*psi2d;
an = (mr*ct2 - ct) ./ ct4;
bn = (ct2 - mr*ct)./ct5;
% cn = ct3./ct5;
% dn = ct3./ct4;
% 
% ax = abs(x);
% ay = abs(y);



%%%%%%%%%% calculate Q and C %%%%%%%%%%
S_FACTOR = exp(-imag(x));
two_n_plus_1=1+2*(1:nmax)';

QS_FAR=S_FACTOR^2*sum(two_n_plus_1.*(abs(an).^2+abs(bn).^2));
CQE= sum(two_n_plus_1.*(conj(psi).*psid - psi.*conj(psid) + bn.*conj(psid).*xi + conj(bn).*psi.*conj(xid) - an.*conj(psi).*xid - conj(an).*conj(xi).*psid));
CQS_NEAR=sum(two_n_plus_1.*(-abs(an).^2.*xid.*conj(xi) + abs(bn).^2.*xi.*conj(xid)));

Qabs = -2*imag((CQE-CQS_NEAR)/m1)/(x0*x0*real(m1));

if Qabs <0 
    Qabs=0;
end
if imag(x) < 10^-6
    QI = 1.0 + (4.0/3.0) * imag(x);
else
    QI = ( ( imag(x) - 0.5)  * exp(2.0 * imag(x)) + 0.5 ) / ( imag(x) ^ 2);
end

Qsca_far=2*QS_FAR/abs(x)/abs(x);
% Qsca_near = 2*imag(-CQS_NEAR/m1)/(n_medium*x0*x0);

SCALE_FACTOR=exp(2*k_medium*x0);
% Csca_scaled=SCALE_FACTOR*Qsca_far*Area;

% Cabs = Qabs*Area;

Qext = (Qabs+Qsca_far) / QI;
w0 = Qsca_far / (Qabs+Qsca_far);

% Qsca=Qext*w0
% Qabs=Qext-Qsca


% QE_INH = (Qabs+Qsca_near) / QI;
% W_INH = Qsca_near / (Qabs+Qsca_near) ;


%%%%%%%%%% calculate g %%%%%%%%%%

xt = 0.0;
for i=1:nmax-1
    xt = xt + (i*(i+2.0))/(i+1.0) * real(an(i)*conj(an(i+1)) + bn(i)*conj(bn(i+1))) + (2.0*i+1.0)/(i*(i+1.0)) * real(an(i)*conj(bn(i)));
end
g_1 = 2 * xt;

xt = 0.0;
for i=1:nmax
    xt = xt + (2*i+1)*(abs(an(i))^2 + abs(bn(i))^2);
end
g = g_1 / xt;



%%%%%%%%%% calculate s1 and s2 %%%%%%%%%%

theta=linspace(eps,pi,nang)';%don't start with zero to avoid division by zero
u=cos(theta);
p=zeros(length(u),nmax);
t=zeros(length(u),nmax);
p(:,1)=1; 
t(:,1)=u;
p(:,2)=3*u; 
t(:,2)=3*cos(2*acos(u));
for n1=3:nmax
    p(:,n1)=((2*n1-1)./(n1-1).*p(:,n1-1).*u)-(n1./(n1-1).*p(:,n1-2));
    t(:,n1)=n1*u.*p(:,n1) - (n1+1).*p(:,n1-1);
end
n2=(2*n+1)./(n.*(n+1));
pin=n2.*p;
tin=n2.*t;
s1=(an'*pin'+bn'*tin');
s2=(an'*tin'+bn'*pin');

%%%%%%%%%% mueller matrix %%%%%%%%%%%%%%%
coeff=0.5*trapz(theta,(abs(s1).^2 + abs(s2).^2).*sin(theta)');
P11 = (abs(s1).^2 + abs(s2).^2)/coeff;
% mm(:,2) = (abs(s2).^2 - abs(s1).^2)/2;
% ct = s1 .* conj(s2);
% mm(:,3)=real(ct);
% mm(:,4)=imag(ct);
