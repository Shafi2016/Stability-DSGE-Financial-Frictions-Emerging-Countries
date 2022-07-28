%% Endogenous discount factor

var  d, c, h, y, i, k, a, lambda,  tb_y, ca_y, r ,util,beta_fun,eta;  

varexo e;                                    
                                             
parameters  gamma, omega, rho, sigmae, delta, psi_1, alpha, phi, beta, r_w, d_bar;
		 alpha  = 0.32;
		  rho    = 0.86;
          phi  = 0.0359;
          r_w    = 0.025;		
         gamma  = 2;
		 omega  = 1.6;
         delta  = 0.03;
         psi_1  = 0.0635; 
		 sigmae = 0.0201;
         beta   = 1/(1+r_w);
		  h_ss   = ((1-alpha)*(alpha/(r_w+delta))^(alpha/(1-alpha)))^(1/(omega-1)); 
		  k_ss   = h_ss/(((r_w+delta)/alpha)^(1/(1-alpha)));
       	 i_ss   = delta*k_ss;                                                     
		 y_ss   = (k_ss^alpha)*(h_ss^(1-alpha));                                   
		c_ss = (1+r_w)^(1/ psi_1) + h_ss^ omega/ omega -1 ;
		tb_ss  = y_ss-c_ss-i_ss;
        d_bar = tb_ss/r_w;
        d_ss   = d_bar;
      
model;
   
         % Evolution of debt
    d = (1+exp(r(-1)))*d(-1)- exp(y)+exp(c)+exp(i)+(phi/2)*(exp(k)-exp(k(-1)))^2;
         % Production function
    exp(y) = exp(a)*(exp(k(-1))^alpha)*(exp(h)^(1-alpha));
         % Law of motion for capital
    exp(k) = exp(i)+(1-delta)*exp(k(-1)); 
          %Euler equation
    exp(lambda)= beta_fun*(1+exp(r))*exp(lambda(+1)); 
         %Definition marginal utility
    exp(lambda)=(exp(c)-((exp(h)^omega)/omega))^(-gamma)-eta*(-psi_1*(1+exp(c)-omega^(-1)*exp(h)^omega)^(-psi_1-1));  
          %Law of motion Lagrange mulitplier on discount factor equation
     eta=-util(+1)+eta(+1)*beta_fun(+1);
             %%Labor FOC
     ((exp(c)-(exp(h)^omega)/omega)^(-gamma))*(exp(h)^(omega-1)) + 

     eta*(-psi_1*(1+exp(c)-omega^(-1)*exp(h)^omega)^(-psi_1-1)*(-exp(h)^(omega-1))) = exp(lambda)*(1-alpha)*exp(y)/exp(h); 
                %Investment FOC
    exp(lambda)*(1+phi*(exp(k)-exp(k(-1)))) = beta_fun*exp(lambda(+1))*(alpha*exp(y(+1))/exp(k)+1-delta+phi*(exp(k(+1))-exp(k))); 
               %Law of motion for TFP
    a = rho*a(-1)+ sigmae*e; 
             %Definition endogenous discount factor   
    beta_fun =(1+exp(c)-omega^(-1)*exp(h)^omega)^(-psi_1);
               %Definition utility function
   util=(((exp(c)-omega^(-1)*exp(h)^omega)^(1-gamma))-1)/(1-gamma);
                %country interest rate
     exp(r) = r_w;
              %Definition of trade balance to ouput ratio
    tb_y = 1-((exp(c)+exp(i)+(phi/2)*(exp(k)-exp(k(-1)))^2)/exp(y));

    ca_y = (1/exp(y))*(d(-1)-d);       
      
end;


initval;
    r     = log((1-beta)/beta);
    d     = d_ss;
    h     = log(h_ss);
    k     = log(k_ss);
    y     = log(y_ss);
    c     = log(c_ss);
    i     = log(i_ss);
    tb_y    = 1-((exp(c)+exp(i))/exp(y));
    util=(((exp(c)-omega^(-1)*exp(h)^omega)^(1-gamma))-1)/(1-gamma);
    beta_fun =(1+exp(c)-omega^(-1)*exp(h)^omega)^(-psi_1);
    eta=-util/(1-beta_fun);
    lambda=log((exp(c)-((exp(h)^omega)/omega))^(-gamma)-eta*(-psi_1*(1+exp(c)-omega^(-1)*exp(h)^omega)^(-psi_1-1)));
end;


check;
steady; 


shocks;
    var e; stderr 1;
   
end;

stoch_simul(order = 3 , irf = 0, periods = 50000, drop = 1000);
//nlma_theoretical_moments = nlma_th_moments(M_,oo_,options_,var_list_);


//simulations = pruning_abounds(M_,options_, 3 , 'lan_meyer-gohde');
  //simulated_moments(M_,oo_, options_, var_list_, 3,  'den_haan_de_wind')
 //simulated_moments(M_,oo_, options_, var_list_, 3,  'andreasen')
 //simulated_moments(M_,oo_, options_, var_list_, 3,  'lan_meyer-gohde')
 //simulated_moments(M_,oo_, options_, var_list_, 3,  'fernandez-villaverde_et_al');
