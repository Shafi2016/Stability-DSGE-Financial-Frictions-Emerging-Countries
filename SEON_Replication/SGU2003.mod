% Endogenous discount factor///// 3rd order perturbation
@#define model1 = 0
% Endogenous discount factor without internalization
@#define model1a = 0
%debt elastic interest rate premium
@#define model2 = 1
%portfolio holding costs
@#define model3 = 0
%Complete markets
@#define model4 =0
  //M_.params(strmatch('psi_1',M_.param_names,'exact')) // USE TO extract psi_1


var  c h y i k a lambda ${\lambda}$ util;  
varexo e;                                    
                                             
parameters  gamma ${\gamma}$
            omega ${\omega}$
            rho ${\gamma}$
            sigma_tfp ${\sigma_{a}}$
            delta ${\delta}$
            psi_1 ${\psi_1}$
            psi_2 ${\psi_2}$
            alpha ${\alpha}$
            phi ${\phi}$
            psi_3 ${\psi_3}$
            psi_4 ${\psi_4}$
            r_bar ${\bar r}$
            d_bar ${\bar d}$;
            
%Table 1

gamma  = 2.61; %risk aversion

omega  = 1.455; %Frisch-elasticity parameter

psi_1  = 0; %set in steady state %elasticity discount factor w.r.t. to arguments of utility function

alpha  = 0.32; %labor share

phi    = 0.0282; %capital adjustment cost parameter

r_bar    = 0.0206; %world interest rate		

delta  = 0.03; %depreciation rate

rho    = 0.8598; %autocorrelation TFP 

sigma_tfp = 0.0186; %standard deviation TFP



%Table 2

psi_2    = 0.000112;

//_bar  = 58;

d_bar  = 100.53;

  d_bar = 77.83;

//psi_3  = 0.0001931;

psi_3  = 0.000166; 

psi_4  = 0; %set in steady state; parameter complete markets case



@#if model1 == 1
var d tb_y, ca_y, r beta_fun, eta;   

model;
    [name='Eq. (4), Evolution of debt']
    d = (1+exp(r(-1)))*d(-1)- exp(y)+exp(c)+exp(i)+(phi/2)*(exp(k)-exp(k(-1)))^2;
    [name='Eq. (5), Production function']
    exp(y) = exp(a)*(exp(k(-1))^alpha)*(exp(h)^(1-alpha));
    [name='Eq. (6), Law of motion for capital']
    exp(k) = exp(i)+(1-delta)*exp(k(-1)); 

    [name='Eq. (8), Euler equation']
    exp(lambda)= beta_fun*(1+exp(r))*exp(lambda(+1)); 
    [name='Eq. (9), Definition marginal utility']
    exp(lambda)=(exp(c)-((exp(h)^omega)/omega))^(-gamma)-eta*(-psi_1*(1+exp(c)-omega^(-1)*exp(h)^omega)^(-psi_1-1));  
    [name='Eq. (10), Law of motion Lagrange mulitplier on discount factor equation']
    eta=-util(+1)+eta(+1)*beta_fun(+1);
    [name='Eq. (11), Labor FOC']
    ((exp(c)-(exp(h)^omega)/omega)^(-gamma))*(exp(h)^(omega-1)) + 
        eta*(-psi_1*(1+exp(c)-omega^(-1)*exp(h)^omega)^(-psi_1-1)*(-exp(h)^(omega-1))) = exp(lambda)*(1-alpha)*exp(y)/exp(h); 
    [name='Eq. (12), Investment FOC']
    exp(lambda)*(1+phi*(exp(k)-exp(k(-1)))) = beta_fun*exp(lambda(+1))*(alpha*exp(y(+1))/exp(k)+1-delta+phi*(exp(k(+1))-exp(k))); 
    [name='Eq. (14), Law of motion for TFP']
    a = rho*a(-1)+sigma_tfp*e; 

    [name='Definition endogenous discount factor, p. 168']
    beta_fun =(1+exp(c)-omega^(-1)*exp(h)^omega)^(-psi_1);
    [name='Definition felicity function']    
    util=(((exp(c)-omega^(-1)*exp(h)^omega)^(1-gamma))-1)/(1-gamma);
    [name='9. Eq. (23), country interest rate']
    exp(r) = r_bar;

    [name='11. p. 169, Definition of trade balance to ouput ratio']
    tb_y = 1-((exp(c)+exp(i)+(phi/2)*(exp(k)-exp(k(-1)))^2)/exp(y));
    ca_y = (1/exp(y))*(d(-1)-d);                                   
end;

steady_state_model;
    r     = log(r_bar);
    d     = d_bar;
    h     = log(((1-alpha)*(alpha/(r_bar+delta))^(alpha/(1-alpha)))^(1/(omega-1)));
    k     = log(exp(h)/(((r_bar+delta)/alpha)^(1/(1-alpha))));
    y     = log((exp(k)^alpha)*(exp(h)^(1-alpha)));
    i     = log(delta*exp(k));
    c     = log(exp(y)-exp(i)-r_bar*d);
  
   tb_y    = 1-((exp(c)+exp(i))/exp(y));
    util=(((exp(c)-omega^(-1)*exp(h)^omega)^(1-gamma))-1)/(1-gamma);
    psi_1=-log(1/(1+r_bar))/(log((1+exp(c)-omega^(-1)*exp(h)^omega)));
    beta_fun =(1+exp(c)-omega^(-1)*exp(h)^omega)^(-psi_1);
    eta=-util/(1-beta_fun);
    lambda=log((exp(c)-((exp(h)^omega)/omega))^(-gamma)-eta*(-psi_1*(1+exp(c)-omega^(-1)*exp(h)^omega)^(-psi_1-1)));
    a     = 0;
    ca_y    = 0;
end;

@# endif

@#if model1a == 1
var d tb_y, ca_y, r beta_fun;   

model;
    [name='Eq. (4), Evolution of debt']
    d = (1+exp(r(-1)))*d(-1)- exp(y)+exp(c)+exp(i)+(phi/2)*(exp(k)-exp(k(-1)))^2;
    [name='Eq. (5), Production function']
    exp(y) = exp(a)*(exp(k(-1))^alpha)*(exp(h)^(1-alpha));
    [name='Eq. (6), Law of motion for capital']
    exp(k) = exp(i)+(1-delta)*exp(k(-1)); 

    [name='Eq. (17), Euler equation']
    exp(lambda)= beta_fun*(1+exp(r))*exp(lambda(+1)); 
    [name='Eq. (18), Definition marginal utility']
    exp(lambda)=(exp(c)-((exp(h)^omega)/omega))^(-gamma);  
    [name='Eq. (11), Labor FOC']
    ((exp(c)-(exp(h)^omega)/omega)^(-gamma))*(exp(h)^(omega-1))= exp(lambda)*(1-alpha)*exp(y)/exp(h); 
    [name='Eq. (12), Investment FOC']
    exp(lambda)*(1+phi*(exp(k)-exp(k(-1)))) = beta_fun*exp(lambda(+1))*(alpha*exp(y(+1))/exp(k)+1-delta+phi*(exp(k(+1))-exp(k))); 
    [name='Eq. (14), Law of motion for TFP']
    a = rho*a(-1)+sigma_tfp*e; 

    [name='Definition endogenous discount factor, p. 168']
    beta_fun =(1+exp(c)-omega^(-1)*exp(h)^omega)^(-psi_1);
    [name='Definition felicity function']    
    util=(((exp(c)-omega^(-1)*exp(h)^omega)^(1-gamma))-1)/(1-gamma);
    [name='Eq. (23), country interest rate']
    exp(r) = r_bar;

    [name='p. 169, Definition of trade balance to ouput ratio']
    tb_y = 1-((exp(c)+exp(i)+(phi/2)*(exp(k)-exp(k(-1)))^2)/exp(y));
    [name='Definition of current account to ouput ratio']
    ca_y = (1/exp(y))*(d(-1)-d);                                   
end;

steady_state_model;
    r     = log(r_bar);
        
    d     = d_bar;
    h     = log(((1-alpha)*(alpha/(r_bar+delta))^(alpha/(1-alpha)))^(1/(omega-1)));
    k     = log(exp(h)/(((r_bar+delta)/alpha)^(1/(1-alpha))));
    y     = log((exp(k)^alpha)*(exp(h)^(1-alpha)));
    i     = log(delta*exp(k));
    c     = log(exp(y)-exp(i)-r_bar*d);
    tb_y    = 1-((exp(c)+exp(i))/exp(y));
    util=(((exp(c)-omega^(-1)*exp(h)^omega)^(1-gamma))-1)/(1-gamma);
    psi_1=-log(1/(1+r_bar))/(log((1+exp(c)-omega^(-1)*exp(h)^omega)));
    beta_fun =(1+exp(c)-omega^(-1)*exp(h)^omega)^(-psi_1);
    lambda= log((exp(c)-((exp(h)^omega)/omega))^(-gamma));
    a     = 0;
    ca_y    = 0;
end;

@# endif

@#if model2 == 1
var d tb_y, ca_y, r riskpremium;
parameters beta ${\beta}$;

model;
    [name='Eq. (4), Evolution of debt']
    d = (1+exp(r(-1)))*d(-1)- exp(y)+exp(c)+exp(i)+(phi/2)*(exp(k)-exp(k(-1)))^2;
    [name='Eq. (5), Production function']
    exp(y) = exp(a)*(exp(k(-1))^alpha)*(exp(h)^(1-alpha));
    [name='Eq. (6), Law of motion for capital']
    exp(k) = exp(i)+(1-delta)*exp(k(-1)); 

    [name='Eq. (24), Euler equation']
    exp(lambda)= beta*(1+exp(r))*exp(lambda(+1)); 
    [name='Eq. (25), Definition marginal utility']
    (exp(c)-((exp(h)^omega)/omega))^(-gamma)   = exp(lambda);  
    [name='Eq. (26), Labor FOC']
    ((exp(c)-((exp(h)^omega)/omega))^(-gamma))*(exp(h)^(omega-1))  = exp(lambda)*(1-alpha)*exp(y)/exp(h); 
    [name='Eq. (27), Investment FOC']
    exp(lambda)*(1+phi*(exp(k)-exp(k(-1)))) = beta*exp(lambda(+1))*(alpha*exp(y(+1))/exp(k)+1-delta+phi*(exp(k(+1))-exp(k))); 
    [name='Eq. (14), Law of motion for TFP']
    a = rho*a(-1)+sigma_tfp*e; 
    [name='Eq. (23), country interest rate']
    exp(r) = r_bar+riskpremium;
    [name='p. 171 below Eq. (28), definition risk premia']
   riskpremium = psi_2*(exp(d-d_bar)-1);
    //riskpremium = psi_2*(exp(d-d_bar)-1)+ (exp(d-d_bar))*d*psi_2; // IDEIR additive  

    [name='p. 169, Definition of trade balance to ouput ratio']
    tb_y = 1-((exp(c)+exp(i)+(phi/2)*(exp(k)-exp(k(-1)))^2)/exp(y));
    ca_y = (1/exp(y))*(d(-1)-d);                                   
    [name='Definition felicity function']    
    util=(((exp(c)-omega^(-1)*exp(h)^omega)^(1-gamma))-1)/(1-gamma);
end;

steady_state_model;
    beta  = 1/(1+r_bar);
    r     = log((1-beta)/beta);
    d     = d_bar;
    h     = log(((1-alpha)*(alpha/(r_bar+delta))^(alpha/(1-alpha)))^(1/(omega-1)));
    k     = log(exp(h)/(((r_bar+delta)/alpha)^(1/(1-alpha))));
    y     = log((exp(k)^alpha)*(exp(h)^(1-alpha)));
    i     = log(delta*exp(k));
    c     = log(exp(y)-exp(i)-r_bar*d);
    tb_y    = 1-((exp(c)+exp(i))/exp(y));
    util=(((exp(c)-omega^(-1)*exp(h)^omega)^(1-gamma))-1)/(1-gamma);
    lambda= log((exp(c)-((exp(h)^omega)/omega))^(-gamma));
    a     = 0;
    ca_y    = 0;
    riskpremium = 0;
end;

@# endif



@#if model3 == 1
var d tb_y, ca_y, r;
parameters beta ${\beta}$;

model;
    [name='Eq. (29), Evolution of debt']
    d = (1+exp(r(-1)))*d(-1)- exp(y)+exp(c)+exp(i)+(phi/2)*(exp(k)-exp(k(-1)))^2+psi_3*(d-d_bar)^2;
    [name='Eq. (5), Production function']
    exp(y) = exp(a)*(exp(k(-1))^alpha)*(exp(h)^(1-alpha));
    [name='Eq. (6), Law of motion for capital']
    exp(k) = exp(i)+(1-delta)*exp(k(-1)); 

    [name='Eq. (30),  Euler equation']
    exp(lambda)*(1-psi_3*(d-d_bar))= beta*(1+exp(r))*exp(lambda(+1)); 
    [name='Eq. (25), Definition marginal utility']
    (exp(c)-((exp(h)^omega)/omega))^(-gamma)   = exp(lambda);  
    [name='Eq. (26), Labor FOC']
    ((exp(c)-((exp(h)^omega)/omega))^(-gamma))*(exp(h)^(omega-1))  = exp(lambda)*(1-alpha)*exp(y)/exp(h); 
    [name='Eq. (27), Investment FOC']
    exp(lambda)*(1+phi*(exp(k)-exp(k(-1)))) = beta*exp(lambda(+1))*(alpha*exp(y(+1))/exp(k)+1-delta+phi*(exp(k(+1))-exp(k))); 
    [name='Eq. (14), Law of motion for TFP']
    a = rho*a(-1)+sigma_tfp*e; 
    [name='Eq. (13), country interest rate']
    exp(r) = r_bar;

    [name='p. 169, Definition of trade balance to ouput ratio']
    tb_y = 1-((exp(c)+exp(i)+(phi/2)*(exp(k)-exp(k(-1)))^2)/exp(y));
    ca_y = (1/exp(y))*(d(-1)-d); 
    [name='Definition felicity function']    
    util=(((exp(c)-omega^(-1)*exp(h)^omega)^(1-gamma))-1)/(1-gamma);
end;

steady_state_model;
    beta  = 1/(1+r_bar);
    r     = log((1-beta)/beta);
    d     = d_bar;
    h     = log(((1-alpha)*(alpha/(r_bar+delta))^(alpha/(1-alpha)))^(1/(omega-1)));
    k     = log(exp(h)/(((r_bar+delta)/alpha)^(1/(1-alpha))));
    y     = log((exp(k)^alpha)*(exp(h)^(1-alpha)));
    i     = log(delta*exp(k));
    c     = log(exp(y)-exp(i)-r_bar*d);
    tb_y    = 1-((exp(c)+exp(i))/exp(y));
    util=(((exp(c)-omega^(-1)*exp(h)^omega)^(1-gamma))-1)/(1-gamma);
    lambda= log((exp(c)-((exp(h)^omega)/omega))^(-gamma));
    a     = 0;
    ca_y    = 0;
end;

@# endif

@#if model4 == 1
    var tb_y;
parameters beta ${\beta}$;

model;
    [name='Eq. (5), Production function']
    exp(y) = exp(a)*(exp(k(-1))^alpha)*(exp(h)^(1-alpha));
    [name='Eq. (6), Law of motion for capital']
    exp(k) = exp(i)+(1-delta)*exp(k(-1)); 
    [name='Eq. (25), Definition marginal utility']
    (exp(c)-((exp(h)^omega)/omega))^(-gamma)   = exp(lambda);  
    [name='Eq. (26), Labor FOC']
    ((exp(c)-((exp(h)^omega)/omega))^(-gamma))*(exp(h)^(omega-1))  = exp(lambda)*(1-alpha)*exp(y)/exp(h); 
    [name='Eq. (27), Investment FOC']
    exp(lambda)*(1+phi*(exp(k)-exp(k(-1)))) = beta*exp(lambda(+1))*(alpha*exp(y(+1))/exp(k)+1-delta+phi*(exp(k(+1))-exp(k))); 
    [name='Eq. (35),  Euler equation']
    exp(lambda)= psi_4; 
    [name='Eq. (14), Law of motion for TFP']
    a = rho*a(-1)+sigma_tfp*e; 
    [name='p. 169, Definition of trade balance to ouput ratio']
    tb_y = 1-((exp(c)+exp(i)+(phi/2)*(exp(k)-exp(k(-1)))^2)/exp(y));
    [name='Definition felicity function']
    util=(((exp(c)-omega^(-1)*exp(h)^omega)^(1-gamma))-1)/(1-gamma);
end;

steady_state_model;
    beta  = 1/(1+r_bar);
    h     = log(((1-alpha)*(alpha/(r_bar+delta))^(alpha/(1-alpha)))^(1/(omega-1)));
    k     = log(exp(h)/(((r_bar+delta)/alpha)^(1/(1-alpha))));
    y     = log((exp(k)^alpha)*(exp(h)^(1-alpha)));
    i     = log(delta*exp(k));
    c     = 0.110602; %from incomplete markets case
    lambda= log((exp(c)-((exp(h)^omega)/omega))^(-gamma));
    psi_4=exp(lambda);
    tb_y    = 1-((exp(c)+exp(i))/exp(y));
    util=(((exp(c)-omega^(-1)*exp(h)^omega)^(1-gamma))-1)/(1-gamma);
    a     = 0;
end;

@# endif


//resid(1);

check;
steady; 


shocks;
   var e; stderr 1;
end;
 stoch_simul(order=3,pruning, periods= 500000,irf= 50);

y_pos=strmatch('y',M_.endo_names,'exact');

c_pos=strmatch('c',M_.endo_names,'exact');

i_pos=strmatch('i',M_.endo_names,'exact');

h_pos=strmatch('h',M_.endo_names,'exact');

tb_y_pos=strmatch('tb_y',M_.endo_names,'exact');

ca_y_pos=strmatch('ca_y',M_.endo_names,'exact');

Z_pos=strmatch('a',M_.endo_names,'exact');
d_pos=strmatch('d',M_.endo_names,'exact');
lambda_pos=strmatch('lambda',M_.endo_names,'exact');

IRF_periods=50;

burnin=5000; %periods for convergence




shock_mat_with_zeros=zeros(burnin+IRF_periods,M_.exo_nbr); %shocks set to 0 to simulate without uncertainty

IRF_no_shock_mat = simult_(oo_.dr.ys,oo_.dr,shock_mat_with_zeros,options_.order)'; %simulate series

stochastic_steady_state=IRF_no_shock_mat(1+burnin,:); % stochastic_steady_state/EMAS is any of the final points after burnin



shock_mat = zeros(burnin+IRF_periods,M_.exo_nbr);

shock_mat(1+burnin,strmatch('e',M_.exo_names,'exact'))= 1;

IRF_mat = simult_(oo_.dr.ys,oo_.dr,shock_mat,options_.order)';



IRF_mat_percent_from_SSS = (IRF_mat(1+burnin+1:1+burnin+IRF_periods,:)-IRF_no_shock_mat(1+burnin+1:1+burnin+IRF_periods,:))./repmat(stochastic_steady_state,IRF_periods,1); %only valid for variables not yet logged



%scale IRFs as reqired

y_vola_IRF 		= 100*IRF_mat_percent_from_SSS(:,y_pos);

c_vola_IRF 		= 100*IRF_mat_percent_from_SSS(:,c_pos);

inv_vola_IRF 	= 100*IRF_mat_percent_from_SSS(:,i_pos);

h_vola_IRF 	= 100*IRF_mat_percent_from_SSS(:,h_pos);

tb_y_vola_IRF 	= 100*IRF_mat_percent_from_SSS(:,tb_y_pos);

ca_y_vola_IRF 	= 100*IRF_mat_percent_from_SSS(:,ca_y_pos);

z_vola_IRF 	= 100*IRF_mat_percent_from_SSS(:,Z_pos);
d_vola_IRF 	= 100*IRF_mat_percent_from_SSS(:,d_pos);
lambda_vola_IRF 	= 100*IRF_mat_percent_from_SSS(:,lambda_pos);

hh=figure;

figure(hh)   

subplot(3,3,1)

hold on

plot(y_vola_IRF,'b-','LineWidth',3)

plot(zeros(IRF_periods,1),'k--','HandleVisibility','off'); xlim([1 IRF_periods]);set(gca,'XTick',[4:4:IRF_periods],'FontSize',12);

title('Output','FontSize',14)

ylabel('Percent','FontSize',12)



figure(hh)   

subplot(3,3,2)

hold on

plot(c_vola_IRF,'b-','LineWidth',3)

plot(zeros(IRF_periods,1),'k--','HandleVisibility','off'); xlim([1 IRF_periods]);set(gca,'XTick',[4:4:IRF_periods],'FontSize',12);

title('Consumption','FontSize',14)

ylabel('Percent','FontSize',12)

%ylim([-0.3 0.1]);set(gca,'YTick',[-0.3:0.1:0.1],'FontSize',12);





figure(hh)   

subplot(3,3,3)

hold on

plot(inv_vola_IRF,'b-','LineWidth',3)

plot(zeros(IRF_periods,1),'k--','HandleVisibility','off'); xlim([1 IRF_periods]);set(gca,'XTick',[4:4:IRF_periods],'FontSize',12);

title('Investment','FontSize',14)

ylabel('Percent','FontSize',12)

%ylim([-0.6 0.4]);set(gca,'YTick',[-0.6:0.2:0.4],'FontSize',12);



figure(hh)   

subplot(3,3,4)

hold on

plot(h_vola_IRF,'b-','LineWidth',3)

plot(zeros(IRF_periods,1),'k--','HandleVisibility','off'); xlim([1 IRF_periods]);set(gca,'XTick',[4:4:IRF_periods],'FontSize',12);

title('Hours worked','FontSize',14)

ylabel('Percent','FontSize',12)



figure(hh)   

subplot(3,3,5)

hold on

plot(tb_y_vola_IRF ,'b-','LineWidth',3)

plot(zeros(IRF_periods,1),'k--','HandleVisibility','off'); xlim([1 IRF_periods]);set(gca,'XTick',[4:4:IRF_periods],'FontSize',12);

title('Trade balance to GDP','FontSize',14)

ylabel('Percent','FontSize',12)



figure(hh)   

subplot(3,3,6)

hold on

plot(ca_y_vola_IRF,'b-','LineWidth',3)

plot(zeros(IRF_periods,1),'k--','HandleVisibility','off'); xlim([1 IRF_periods]);set(gca,'XTick',[4:4:IRF_periods],'FontSize',12);

title('current accout to GDP','FontSize',14)

ylabel('Percent','FontSize',12)

figure(hh)   

subplot(3,3,7)

hold on

plot(z_vola_IRF,'b-','LineWidth',3)

plot(zeros(IRF_periods,1),'k--','HandleVisibility','off'); xlim([1 IRF_periods]);set(gca,'XTick',[4:4:IRF_periods],'FontSize',12);

title('z','FontSize',14)

ylabel('Percent','FontSize',12)

subplot(3,3,8)

hold on

plot(d_vola_IRF ,'b-','LineWidth',3)

plot(zeros(IRF_periods,1),'k--','HandleVisibility','off'); xlim([1 IRF_periods]);set(gca,'XTick',[4:4:IRF_periods],'FontSize',12);

title('d','FontSize',14)

ylabel('Percent','FontSize',12)

subplot(3,3,9)

hold on

plot(lambda_vola_IRF,'b-','LineWidth',3)

plot(zeros(IRF_periods,1),'k--','HandleVisibility','off'); xlim([1 IRF_periods]);set(gca,'XTick',[4:4:IRF_periods],'FontSize',12);

title('Lambda','FontSize',14)

ylabel('Percent','FontSize',12)


