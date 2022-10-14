var Y I K N A C W r m pi i P M ;

varexo e em;

parameters alpha beta zeta psi delta theta rho rhom CONS sigmae pistar;
alpha = 1/3;
beta = 0.99;    
zeta=2;
psi=1;
delta = 0.025;
theta = 4;
rho = 0.97;
rhom=0.5;
CONS = (alpha/(1/beta-(1-delta)))^(1/(1-alpha));
sigmae = 0.01;
pistar=0;


model;
theta /(1-N) = 1/C*W;
1/(C)=beta/C(+1)*(1+r);
(r) = alpha*(A)*(K(-1))^(alpha-1)*N^(1-alpha)-delta;
(W) = (1-alpha)*(A)*(K(-1))^(alpha)*N^(-alpha);
(Y) = (A)*(K(-1))^(alpha)*N^(1-alpha);
K = I + (1-delta)*K(-1);
(Y) = (C) + (I);
ln(m) - ln(m(-1)) = (1-rhom)*pistar -pi+ rhom*pi(-1) + rhom*(ln(m(-1)) - ln(m(-2))) + em;
m= psi^zeta * C^zeta * (1+1/i) ^zeta;
1+r = (1+i)/(1+pi(+1)); 
ln(A) = rho*ln(A(-1))+ e ;
ln(P) = pi + ln(P(-1));
ln(M)=ln(m)*ln(P);

end;

initval;
N= (1-alpha)/theta * CONS^alpha / ( (theta+1-alpha)/theta * CONS^alpha -delta * CONS);K = CONS *N;
K = CONS *N;
Y = (K/N)^alpha*N;
C = N* ((K/N)^alpha - delta*(K/N));
I = Y-C;
A = 1;
r = 1/beta-1;
W = (1-alpha)*(K/N)^alpha;
i=1/beta-1;
pi=0;
P=0.001;
m= psi^zeta * C^zeta * (1+1/i) ^zeta;
M=m*P;
end;


shocks;
var e = sigmae^2;
var em = sigmae^2;
end;

STEADY;

set_dynare_seed=7;
stoch_simul(order=1 , irf=40);
