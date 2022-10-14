var Y I K N A C G W r;
varexo ea eg;

parameters alpha beta chi delta theta rho rhog omega CONS sigmae;
alpha = 1/3;
beta = 0.99;
chi=1;
delta = 0.025;
theta = 4; 
rho = 0.97; 
rhog=0.95; 
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
(r) = alpha* (A) * (K(-1))^(alpha-1) *N*(1-alpha)-delta;
(W) = (1-alpha) * (A)*(K(-1)) ^ (alpha) *N^(-alpha);

end;

x_array = [0.95,0.7,0.99];
for count = 1:3
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
 end;


STEADY;
set_dynare_seed=7;
stoch_simul(order=2, irf=40);

Yi(:,count)=Y_eg;
Ii(:,count)=I_eg;
Ki(:,count)=K_eg;
Ni(:,count)=N_eg;
Ai(:,count)=A_eg;
Ci(:,count)=C_eg;
Wi(:,count)=W_eg;
ri(:,count)=r_eg;
Gi(:,count)=G_eg;

end;


t = linspace(0,40,40);
for i =1:40
line(i) = 0;
end
line=line.';
figure
subplot(3,3,1)
plot(t,line,'red',t,Yi(:,2),'black',t,Yi(:,3),'--red',t,Yi(:,1),'--blue','LineWidth',1)
title("(Y)")
xlabel('t')
ylabel('value')
legend('Y=0','rho=0.95','rho=0.99','rho=0.5')
subplot(3,3,2)
plot(t,line,'red',t,Ii(:,2),'black',t,Ii(:,3),'--red',t,Ii(:,1),'--blue','LineWidth',1)
title("(I)")
xlabel('t')
ylabel('value')
subplot(3,3,3)
plot(t,line,'red',t,Ki(:,2),'black',t,Ki(:,3),'--red',t,Ki(:,1),'--blue','LineWidth',1)
title("(K)")
xlabel('t')
ylabel('value')
subplot(3,3,4)
plot(t,line,'red',t,Ni(:,2),'black',t,Ni(:,3),'--red',t,Ni(:,1),'--blue','LineWidth',1)
title("(N)")
xlabel('t')
ylabel('value')
subplot(3,3,5)
plot(t,line,'red',t,Gi(:,2),'black',t,Gi(:,3),'--red',t,Gi(:,1),'--blue','LineWidth',1)
title("(G)")
xlabel('t')
ylabel('value')
subplot(3,3,6)
plot(t,line,'red',t,Ci(:,2),'black',t,Ci(:,3),'--red',t,Ci(:,1),'--blue','LineWidth',1)
title("C")
xlabel('t')
ylabel('value')
subplot(3,3,7)
plot(t,line,'red',t,Wi(:,2),'black',t,Wi(:,3),'--red',t,Wi(:,1),'--blue','LineWidth',1)
title("(W)")
xlabel('t')
ylabel('value')
subplot(3,3,8)
plot(t,line,'red',t,ri(:,2),'black',t,ri(:,3),'--red',t,ri(:,1),'--blue','LineWidth',1)
title("r")
xlabel('t')
ylabel('value')



 



