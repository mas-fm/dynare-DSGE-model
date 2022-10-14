var Y I K N A C G W r Tk Tn;
varexo ea eg en ek;

parameters alpha beta chi delta theta rho rhog rhok rhon tk tn omega CONS sigmae;
alpha = 1/3;
beta = 0.99;
chi=1;
delta = 0.025;
theta = 4; 
rho = 0.97; 
rhog=0.95; 
rhok=0.9; 
rhon=0.9; 
tn=0.2;
tk=0.1;
omega = 0.2;
CONS = (alpha/(1/beta-(1-delta)))^(1/(1-alpha)); //CONS=k/N IN STEADY STATE
sigmae = 0.01;

model;
theta * N^chi = 1/C*W;
1/(C)=beta / (C(+1))*(1+r);
(K) = I + (1-delta) * (K(-1));
ln(A) = rho*ln(A(-1)) + ea;
ln(G) = (1-rhog) * ln (omega * Y) + rhog *ln (G(-1)) +eg;
(Y) = (A)*(K(-1))^(alpha) *N*(1-alpha);
(Y) = (C) + (I) + G;
(r) = (1-Tk) *alpha* (A)*(K(-1))^(alpha-1) *N*(1-alpha)-delta;
(W) = (1-Tn)*(1-alpha) * (A)*(K(-1))^(alpha) *N^(-alpha);
Tk=(1-rhok)*tk+rhok* Tk(-1) +ek;
Tn=(1-rhon) *tn+rhon* Tn(-1) +en;
end;

x_array = [0.97];
for count = 1:1
    rhog=x_array(count);

initval;
 N=(1/theta * ((1-alpha) * (CONS)^alpha) / ((CONS)^alpha* (1-omega) -delta* (CONS)) )^(1/(1+chi));
 K = CONS *N;
 Y = (K/N)^alpha*N;
 C = N* ((1-omega) * (K/N)^alpha - delta* (K/N));
 I = Y-C-G;
 A = 1;
 G= omega *Y;
 r = 1/beta-1;
 W = (1-alpha) * (K/N)^alpha;
 end;
  
 shocks;
 var ea = sigmae^2;
 var eg = sigmae^2;
 var en = sigmae^2;
 var ek = sigmae^2;
 end;


STEADY;
stoch_simul(order=2, irf=100);

end;
