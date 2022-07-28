% Portfolio adjustment cost model 

var  d, c, h, y, i, k, a, lambda,  tb, ca, r ;  

varexo e;                                    
                                             
parameters  gamma, omega, rho, sigmae, delta, psi, alpha, phi, beta, r_w, d_bar;
		 alpha  = 0.32;
		  rho    = 0.86;
          phi  = 0.0359;
          r_w    = 0.025;		
         gamma  = 2;
		 omega  = 1.6;
         delta  = 0.03;
            //psi = 0.00006;// lan_meyer-gohde
            //psi = 0.00015;//fernandez-villaverde_et_al
            // psi = 0.00008;//andreasen pruning
          // psi    = 0.00009;// Den Haan Pruning 
           //psi    = 0.00008;//with out pruning
              psi = 0.00007; //NLMA
         sigmae = 0.0201;
         beta   = 1/(1+r_w);
		  h_ss   = ((1-alpha)*(alpha/(r_w+delta))^(alpha/(1-alpha)))^(1/(omega-1)); 
		  k_ss   = h_ss/(((r_w+delta)/alpha)^(1/(1-alpha)));
       	  i_ss   = delta*k_ss;                                                     
		 y_ss   = (k_ss^alpha)*(h_ss^(1-alpha));                                   
		  d_bar  = 57.75;
          d_ss   = d_bar;                                                        
		c_ss   = y_ss-i_ss-r_w*d_ss;
		tb_ss  = y_ss-c_ss-i_ss;
        
       
model;
    d = (1+exp(r(-1)))*d(-1)- exp(y)+exp(c)+exp(i)+(phi/2)*(exp(k)-exp(k(-1)))^2+psi*(d-d_bar)^2;
    exp(y) = exp(a)*(exp(k(-1))^alpha)*(exp(h)^(1-alpha));
    exp(k) = exp(i)+(1-delta)*exp(k(-1));    
    exp(lambda)*(1-psi*(d-d_bar))= beta*(1+exp(r))*exp(lambda(+1)); 
    (exp(c)-((exp(h)^omega)/omega))^(-gamma)   = exp(lambda);  
  ((exp(c)-((exp(h)^omega)/omega))^(-gamma))*(exp(h)^omega)  = exp(lambda)*(1-alpha)*exp(y); 
    exp(lambda)*(1+phi*(exp(k)-exp(k(-1)))) = beta*exp(lambda(+1))*(alpha*exp(y(+1))/exp(k)+1-delta+phi*(exp(i(+1))-delta*exp(k))); 
    a = rho*a(-1)+ sigmae*e;
    tb = 1-((exp(c)+exp(i))/exp(y));
    ca = (1/exp(y))*(d-d(-1));                                   
    exp(r) = r_w; 
      
end;


initval;
    r     = log((1-beta)/beta);
    d     = d_ss;
    h     = log(h_ss);
    k     = log(k_ss);
    y     = log(y_ss);
    c     = log(c_ss);
    i     = log(i_ss);
    tb    = 1-((exp(c)+exp(i))/exp(y));
    lambda= log((exp(c)-((exp(h)^omega)/omega))^(-gamma));
end;


check;
steady; 


shocks;
    
   
    var e; stderr 1;
end;



stoch_simul(order = 3 , irf = 0 , periods = 0, drop = 1000);
nlma_theoretical_moments = nlma_th_moments(M_,oo_,options_,var_list_);
//simulations = pruning_abounds(M_,options_, 3 , 'lan_meyer-gohde');
  //simulated_moments(M_,oo_, options_, var_list_, 3,  'den_haan_de_wind');
 //simulated_moments(M_,oo_, options_, var_list_, 3,  'andreasen');
  //simulated_moments(M_,oo_, options_, var_list_, 3,  'lan_meyer-gohde');
 //simulated_moments(M_,oo_, options_, var_list_, 3,  'fernandez-villaverde_et_al');