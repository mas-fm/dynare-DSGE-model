var Y I K N A C W r R m pi i P lambda mu dm;

varexo em;

parameters alpha beta zeta psi delta theta rho rhom CONS sigmae pistar;
alpha = 1/3;
beta = 0.99;
zeta=2;
psi=1;
delta = 0.025;
theta = 4;
rho = 0.97;
rhom=0.5;
CONS = (alpha/(1/beta-(1-delta)))^(1/(1-alpha)); //CONS=k/N IN STEADY STATE
sigmae = 0.01;
pistar=0;


model;
theta /(1-N) = lambda*W;//1
1/(C)=lambda+mu;//2
lambda= beta * (lambda(1)*(1+i)*P/P(+1) );//
lambda= beta * ( lambda(+1)*(1+R(1)-delta)  );//
lambda= beta * (mu(1) +lambda(1))*P/P(1);//
(K) = I + (1-delta)*(K(-1));//6
ln(A) = rho*ln(A(-1)) ;//7
(Y) = (A)*(K(-1))^(alpha)*N^(1-alpha);//
(Y) = (C) + (I);//9
(R) = alpha*(A)*(K(-1))^(alpha-1)*N^(1-alpha);//
(W) = (1-alpha)*(A)*(K(-1))^(alpha)*N^(-alpha);//
dm = (1-rhom)*pistar -pi + rhom*pi(-1) + rhom*dm(-1)+em;//
1+r = (1+i)/(1+pi(+1))^(-1);  //13
ln(P) = pi + ln(P(-1));//14
//ln(M) = ln(m) + ln(P);
m=C;//
dm = ln(m) - ln(m(-1));
end;

////////////////////////////
/////////////////////////////////////////////
//x_array = [0.97];
//for count = 1:3
//    rhog=x_array(count);
/////////////////////////////////////////////
////////////////////////////


initval;
    N = ( 1+ (1-delta*CONS^(1-alpha)) *theta/beta/(1-alpha) ) ^(-1);
//N=2.5;
K = CONS *N;
Y = (K/N)^alpha*N;
C = N* ((K/N)^alpha - delta*(K/N));
I = Y-C;
A = 1;
r = 1/beta-1;
R=alpha * CONS^(alpha-1);
W = (1-alpha)*(K/N)^alpha;
i=1/beta-1;
pi=0;
P=1;
lambda=beta/C;
mu=lambda*(1/beta -1);
m=C;
//M=P*m;
end;


shocks;
var em = sigmae^2;
end;

STEADY;



//Solve for the stochastic dynamics, 7 as for ps7
set_dynare_seed=7;
stoch_simul;


/*
y_s(:,count)=Y_eg;
i_s(:,count)=I_eg;
k_s(:,count)=K_eg;
n_s(:,count)=N_eg;
a_s(:,count)=G_eg;
c_s(:,count)=C_eg;
w_s(:,count)=W_eg;
r_s(:,count)=r_eg;


end;


// Display the path of consumption and capital
qq=40;
t = linspace(0,qq,qq);
t=t.';
for i =1:qq
line(i) = 0;
end
line=line.';


figure
subplot(3,3,1)
plot(t,line,'red',t,y_s(:,2),'black',t,y_s(:,3),'--red',t,y_s(:,1),'--blue','LineWidth',1)
title("Y")
xlabel('t')
ylabel('value')
legend('y=0','rho=0.7','rho=0.99','rho=0.95')

subplot(3,3,2)
plot(t,line,'red',t,k_s(:,2),'black',t,k_s(:,3),'--red',t,k_s(:,1),'--blue','LineWidth',1)
title("K")
xlabel('t')
ylabel('value')

subplot(3,3,3)
plot(t,line,'red',t,c_s(:,2),'black',t,c_s(:,3),'--red',t,c_s(:,1),'--blue','LineWidth',1)
title("C")
xlabel('t')
ylabel('value')

subplot(3,3,4)
plot(t,line,'red',t,n_s(:,2),'black',t,n_s(:,3),'--red',t,n_s(:,1),'--blue','LineWidth',1)
title("N")
xlabel('t')
ylabel('value')

subplot(3,3,5)
plot(t,line,'red',t,w_s(:,2),'black',t,w_s(:,3),'--red',t,w_s(:,1),'--blue','LineWidth',1)
title("W")
xlabel('t')
ylabel('value')

subplot(3,3,6)
plot(t,line,'red',t,r_s(:,2),'black',t,r_s(:,3),'--red',t,r_s(:,1),'--blue','LineWidth',1)
title("r")
xlabel('t')
ylabel('value')

subplot(3,3,7)
plot(t,line,'red',t,i_s(:,2),'black',t,i_s(:,3),'--red',t,i_s(:,1),'--blue','LineWidth',1)
title("I")
xlabel('t')
ylabel('value')

subplot(3,3,8)
plot(t,line,'red',t,a_s(:,2),'black',t,a_s(:,3),'--red',t,a_s(:,1),'--blue','LineWidth',1)
title("G")
xlabel('t')
ylabel('value')
*/


