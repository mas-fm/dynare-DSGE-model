var n w c m i pi mc a y vp x1 x2 pisharp x yf;

varexo ea ei;

parameters beta sigmae epsilon pistar istar phi phipi phix theta rhoa rhoi chi xstar ;

beta =0.99;
sigmae=0.01;
epsilon=10;
pistar=0;
istar=(1/beta)-1;
phi=0.75;
phipi=1.5;
phix=0;
theta=1;
rhoa=0.95;
rhoi=0.8;
chi=0.5;
xstar=1;


model;

theta/(1-n) = w/c;
m^(-chi) = (1/c)*(i/(1+i));
1/c = beta*(1+i)/(c(+1)*(1+pi(+1)));
mc = w/a;
c=y;
y=a*n/vp;
vp=(1-phi)*(1+pisharp)^(-epsilon)*(1+pi)^(epsilon)+(1+pi)^(epsilon)*phi*vp(-1);
(1+pi)^(1-epsilon) = (1-phi)*((1+pisharp)^(1-epsilon))+phi ;
x1 = c^(-1) *mc*y+ phi *beta*(1+pi(1))^(epsilon)*x1(1) ;
x2 = c^(-1)*   y + phi *beta*(1+pi(1))^(epsilon-1)*x2(1)  ;
1+pisharp = (epsilon/(epsilon-1)) * (1+pi) *x1/x2;
ln(a) = rhoa*ln(a(-1))+ea;
ln(x)=ln(n/vp)-ln((epsilon-1)/(epsilon-1+(epsilon*theta)));
i = (1-rhoi)*istar + rhoi*(i(-1))+(1-rhoi)*(phipi*(pi-pistar) + phix*(ln(x)-ln(xstar))) + ei ;
yf = a*((1+(theta*epsilon)/(epsilon-1)))^(-1);
end;
/////////
phi_array = linspace(0.1,0.99,10);
for count =1:10
    phi=phi_array(count);

////
initval;
a=1;

i=1/beta*(1+pistar)-1;
vp = 1;
x1 = mc*x2;
x2 =3.88;
mc = (epsilon-1)/epsilon;
w=mc;
n = (mc*vp)/(mc*vp+theta);
y=n/vp;
m = ((1/y)*(i/(1+i)))^(-1/chi);
c=y;
x=xstar;
pi=0;
pisharp=0;

end;

shocks;
var ea = sigmae^2;
var ei = sigmae^2;
end;


steady;
n_s(1,count) = round(oo_.steady_state(1,1))
w_s(1,count) =round (oo_.steady_state(2,1))
c_s(1,count) =round (oo_.steady_state(3,1))
m_s(1,count) = round (oo_.steady_state(4,1))
i_s(1,count) =round (oo_.steady_state(5,1))
pi_s(1,count) = round (oo_.steady_state(6,1))
mc_s(1,count) =round (oo_.steady_state(7,1))
a_s(1,count) =round (oo_.steady_state(8,1))
y_s(1,count) =round (oo_.steady_state(9,1))
vp_s(1,count) =round (oo_.steady_state(10,1))
x1_s(1,count) =round (oo_.steady_state(11,1))
x2_s(1,count) =round( oo_.steady_state(12,1))
pisharp_s(1,count) =round (oo_.steady_state(13,1))
x_s(1,count) =round (oo_.steady_state(14,1))
yf_s(1,count) =round (oo_.steady_state(15,1))
phi_s(1,count) = phi
end;

figure
subplot(3,5,1)
plot(phi_s(1,:),n_s(1,:))
xlabel('phi')
ylabel('value')
title("N")


subplot(3,5,2)
plot(phi_s(1,:),w_s(1,:))
xlabel('phi')
ylabel('value')
title("W")

subplot(3,5,3)
plot(phi_s(1,:),c_s(1,:))
xlabel('phi')
ylabel('value')
title("C")


subplot(3,5,4)
plot(phi_s(1,:),m_s(1,:))
xlabel('phi')
ylabel('value')
title("m")


subplot(3,5,5)
plot(phi_s(1,:),i_s(1,:))
xlabel('phi')
ylabel('value')
title("i")


subplot(3,5,6)
plot(phi_s(1,:),pi_s(1,:))
xlabel('phi')
ylabel('value')
title("pi")


subplot(3,5,7)
plot(phi_s(1,:),mc_s(1,:))
xlabel('phi')
ylabel('value')
title("mc")

subplot(3,5,8)
plot(phi_s(1,:),a_s(1,:))
xlabel('phi')
ylabel('value')
title("a")


subplot(3,5,9)
plot(phi_s(1,:),y_s(1,:))
xlabel('phi')
ylabel('value')
title("Y")


subplot(3,5,10)
plot(phi_s(1,:),vp_s(1,:))
xlabel('phi')
ylabel('value')
title("Vp")


subplot(3,5,11)
plot(phi_s(1,:),x1_s(1,:))
xlabel('phi')
ylabel('value')
title("x1")


subplot(3,5,12)
plot(phi_s(1,:),x2_s(1,:))
xlabel('phi')
ylabel('value')
title("x2")


subplot(3,5,13)
plot(phi_s(1,:),pisharp_s(1,:))
xlabel('phi')
ylabel('value')
title("Psharp")
