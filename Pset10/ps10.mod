var Y N A C W mc v m i pi pisha dm x1 x2;

varexo e em;

parameters  beta eta psi epsilon theta rho rhom sigmae pistar phi sigma;

beta = 0.99;
eta=1;
psi=1;
epsilon = 10;
theta = 1;
rho = 0.95;
rhom=0.0;
sigmae = 0.01;
pistar = 0;
phi=0.75;
sigma=1;

model;
C^(-sigma) = beta* C(1)^(-sigma) * (1+i) /( 1+pi(1) ) ;//1
psi * N^eta = C^(-sigma)*W ;//2
m = theta * (1+i)/i * C^(sigma);//3
mc= W/A;//4
C = Y; //5
Y=A*N/v; //6
v=(1-phi)*(1+pisha)^(-epsilon)*(1+pi)^(epsilon)+(1+pi)^(epsilon)*phi*v(-1); //7
(1+pi)^(1-epsilon) = (1-phi)*(1+pisha)^(1-epsilon)+phi ;//8
1+pisha = epsilon/(epsilon-1) * (1+pi) *x1/x2 ;//9
x1 = C^(-sigma) *mc*Y+ phi *beta*(1+pi(1))^(epsilon)*x1(1) ;//10
x2 = C^(-sigma)*   Y + phi *beta*(1+pi(1))^(epsilon-1)*x2(1)  ;//11
ln(A)= rho*ln(A(-1)) + e; //12
dm = (1-rhom)*pistar -pi + rhom*pi(-1) + rhom*dm(-1)+em;//13
dm = ln(m) - ln(m(-1));
end;



initval;
A = 1;
dm=0;
i=1/beta*(1+pistar)-1;
pisha=0;
pi=0;
v=1;
W=mc;
x2=10;
x1=x2*mc;
mc=(epsilon-1)/epsilon;
N = ( 1/psi* v^sigma * mc )^(1/(eta+sigma));
Y = N/v;
C = Y;
m=theta*(1+i)/i*Y^sigma;

end;

shocks;
var e = sigmae^2;
var em = sigmae^2;
end;


STEADY;



//Solve for the stochastic dynamics, 10 as for ps1
set_dynare_seed=7;
stoch_simul(order=1,irf=100);
