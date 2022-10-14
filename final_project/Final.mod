var pi           
    y_gap        
    y_nat      
    y      
    yhat     
    r_nat       
    r_real         
    i       
    n           
    m_real     
    m_nominal          
    a           
    z           
    p          
    w           
    c          
    w_real  

    money_growth
    mu
    mu_hat
;     

varexo  eps_a
        eps_m  
        eps_z   
       ;

parameters alppha      
    betta              
    rho_a              
    rho_m       
    rho_z   
    siggma   
    varphi                    
    eta        
    epsilon             
    theta  
    ;


siggma = 1;
varphi=5;
theta=3/4;
rho_m =0.5;
rho_z  = 0.5;
rho_a  = 0.9;
betta  = 0.99;
eta  =3.77;
alppha=1/4;
epsilon=9;

%---------------------------------------------------------------
%  
%---------------------------------------------------------------

model; 

#Omega=(1-alppha)/(1-alppha+alppha*epsilon);       
#psi_n_ya=(1+varphi)/(siggma*(1-alppha)+varphi+alppha); 
#lambda=(1-theta)*(1-betta*theta)/theta*Omega;     
#kappa=lambda*(siggma+(varphi+alppha)/(1-alppha));     

[name='New Keynesian Phillips Curve eq. (22)']
pi = betta*pi(+1)+kappa*y_gap;

[name='Dynamic IS Curve eq. (23)']
y_gap=-1/siggma*(i-pi(+1)-r_nat)+y_gap(+1);


[name='Definition natural rate of interest eq. (24)']
r_nat=-siggma*psi_n_ya*(1-rho_a)*a+(1-rho_z)*z;

[name='Definition real interest rate']
r_real=i-pi(+1);

[name='Definition natural output, eq. (20)']
y_nat = psi_n_ya*a;

[name='Definition output gap']
y_gap=y-y_nat;

[name='TFP shock']
a = rho_a*a(-1) + eps_a;

[name='Production function (eq. 14)']
y = a+(1-alppha)*n;

[name='Preference shock, p. 54']
z     = rho_z*z(-1) - eps_z;

[name='Real money demand (eq. 4)']
m_real=y-eta*i;

[name='definition nominal money growth']
money_growth = m_real - m_real(-1)+pi ;

[name='exogenous process for money supply growth rate']
money_growth = rho_m * money_growth(-1) + eps_m;

[name='Output deviation from steady state']
yhat=y-steady_state(y);

[name='Definition price level']
pi = p - p(-1);

[name='resource constraint, eq. (12)']
y=c;

[name='FOC labor, eq. (2)']
w-p=siggma*c+varphi*n;

[name='definition real wage']
w_real= w - p;

[name='definition nominal money stock']
m_nominal = m_real + p;

[name='average price markup, eq. (18)']
mu = -(siggma+(varphi+alppha)/(1-alppha))*y+(1+varphi)/(1-alppha)*a;

[name='average price markuo, eq. (20)']
mu_hat=-(siggma+(varphi+alppha)/(1-alppha))*y_gap;


end;


%----------------------------------------------------------------
%  Part 1 , 2
%---------------------------------------------------------------


shocks;
    var eps_a = 1^2; 
end;

steady;
check;

%---------------------------------------------------------------
%  Part 1 , 2
%---------------------------------------------------------------


stoch_simul(order = 1,irf=15,noprint,nograph);

hFig = figure(1); clf;
set(hFig, 'Position', [0 0 1000 800])
subplot(5,2,1); plot(y_gap_eps_a, '-o'); title('Output gap');
subplot(5,2,2); plot(4*pi_eps_a, '-o'); title('Inflation');
subplot(5,2,3); plot(y_eps_a, '-o'); title('Output');
subplot(5,2,4); plot(n_eps_a, '-o'); title('Employment');
subplot(5,2,5); plot(w_real_eps_a, '-o'); title('Real wage');
subplot(5,2,6); plot(p_eps_a, '-o'); title('Price level');
subplot(5,2,7); plot(4*i_eps_a, '-o'); title('Nominal rate');
subplot(5,2,8); plot(4*r_real_eps_a, '-o'); title('Real rate');
subplot(5,2,9); plot(m_nominal_eps_a, '-o'); title('Money supply');
subplot(5,2,10);plot(a_eps_a, '-o'); title('a');
saveas(figure(1) , 'Figures1.eps' , 'epsc')


%----------------------------------------------------------------
%  Part 3 
%---------------------------------------------------------------

theta = 0.9;

shocks;
    var eps_a = 1^2; 
end;

stoch_simul(order = 1,irf=30,noprint,nograph) ;

t = zeros(10,30);


t(1,:)= y_gap_eps_a;
t(2,:)= 4*pi_eps_a;
t(3,:)= y_eps_a;
t(4,:)= n_eps_a;
t(5,:)= w_real_eps_a;
t(6,:)= p_eps_a;
t(7,:)= 4*i_eps_a;
t(8,:)= 4*r_real_eps_a;
t(9,:)= m_nominal_eps_a;
t(10,:)= a_eps_a;

theta = 0.1;

shocks;
    var eps_a = 1^2; 
end;

stoch_simul(order = 1,irf=30,noprint,nograph) ;

hFig = figure(2); clf;
set(hFig, 'Position', [0 0 1000 1000])
subplot(5,2,1); plot(1:30,y_gap_eps_a, '-o',1:30,t(1,:),'r-o'); title('Output gap');legend('\theta = 0.1' ,'\theta = 0.9');
subplot(5,2,2); plot(1:30,4*pi_eps_a, '-o',1:30,t(2,:),'r-o'); title('Inflation');legend('\theta = 0.1' ,'\theta = 0.9');
subplot(5,2,3); plot(1:30,y_eps_a, '-o',1:30,t(3,:),'r-o'); title('Output');legend('\theta = 0.1' ,'\theta = 0.9');
subplot(5,2,4); plot(1:30,n_eps_a, '-o',1:30,t(4,:),'r-o'); title('Employment');legend('\theta = 0.1' ,'\theta = 0.9');
subplot(5,2,5); plot(1:30,w_real_eps_a, '-o',1:30,t(5,:),'r-o'); title('Real wage');legend('\theta = 0.1' ,'\theta = 0.9');
subplot(5,2,6); plot(1:30,p_eps_a, '-o',1:30,t(6,:),'r-o'); title('Price level');legend('\theta = 0.1' ,'\theta = 0.9');
subplot(5,2,7); plot(1:30,4*i_eps_a, '-o',1:30,t(7,:),'r-o'); title('Nominal rate');legend('\theta = 0.1' ,'\theta = 0.9');
subplot(5,2,8); plot(1:30,4*r_real_eps_a, '-o',1:30,t(8,:),'r-o'); title('Real rate');legend('\theta = 0.1' ,'\theta = 0.9');
subplot(5,2,9); plot(1:30,m_nominal_eps_a, '-o',1:30,t(9,:),'r-o'); title('Money supply');legend('\theta = 0.1' ,'\theta = 0.9');
subplot(5,2,10); plot(1:30,a_eps_a, '-o',1:30,t(10,:),'r-o'); title('a');legend('\theta = 0.1' ,'\theta = 0.9');
saveas(figure(2) , 'Figures2.eps' , 'epsc')



%---------------------------------------------------------------
%  Part 4,5 
%---------------------------------------------------------------


@#define thetanumb = 100

@#define irfnum = 30

thetavals = linspace(0.01,.99,@{thetanumb}) ;

thetanumb  = @{thetanumb};

result = zeros(1,10);

@#for step in 1:thetanumb

theta =  thetavals(@{step});

shocks;
    var eps_a = 1^2; 
end;

stoch_simul(order = 1,irf=@{irfnum},noprint,nograph) ;

t = zeros(9,@{irfnum});

t(1,:)= theta;
t(2,:)= y_gap_eps_a;
t(3,:)= 4*pi_eps_a;
t(4,:)= y_eps_a;
t(5,:)= n_eps_a;
t(6,:)= w_real_eps_a;
t(7,:)= p_eps_a;
t(8,:)= 4*i_eps_a;
t(9,:)= 4*r_real_eps_a;
t(10,:)= m_nominal_eps_a;
t(11,:)= m_nominal_eps_a;
t(12,:)= m_nominal_eps_a;


devi = max(abs(t)');
result(@{step},:) = devi;

@#endfor


hFig = figure(3); clf;
set(hFig, 'Position', [0 0 1000 800])
subplot(5,2,1); plot(result(:,1),result(:,2), '-'); title('Output gap');xlabel('\theta'); 
subplot(5,2,2); plot(result(:,1),result(:,3), '-'); title('Inflation');xlabel('\theta');
subplot(5,2,3); plot(result(:,1),result(:,4), '-'); title('Output');xlabel('\theta');
subplot(5,2,4); plot(result(:,1),result(:,5), '-'); title('Employment');xlabel('\theta');
subplot(5,2,5); plot(result(:,1),result(:,6), '-'); title('Real wage');xlabel('\theta');
subplot(5,2,6); plot(result(:,1),result(:,7), '-'); title('Price level');xlabel('\theta');
subplot(5,2,7); plot(result(:,1),result(:,8), '-'); title('Nominal rate');xlabel('\theta');
subplot(5,2,8); plot(result(:,1),result(:,9), '-'); title('Real rate');xlabel('\theta');
subplot(5,2,9); plot(result(:,1),result(:,10),'-'); title('Money supply');xlabel('\theta');
saveas(figure(3) , 'Figures3.eps' , 'epsc')

