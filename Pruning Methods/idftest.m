%
% Status : main Dynare file 
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

clear all
clear global
tic;
global M_ oo_ options_ ys0_ ex0_ estimation_info
options_ = [];
M_.fname = 'idftest';
%
% Some global variables initialization
%
global_initialization;
diary off;
logname_ = 'idftest.log';
if exist(logname_, 'file')
    delete(logname_)
end
diary(logname_)
M_.exo_names = 'e';
M_.exo_names_tex = 'e';
M_.endo_names = 'd';
M_.endo_names_tex = 'd';
M_.endo_names = char(M_.endo_names, 'c');
M_.endo_names_tex = char(M_.endo_names_tex, 'c');
M_.endo_names = char(M_.endo_names, 'h');
M_.endo_names_tex = char(M_.endo_names_tex, 'h');
M_.endo_names = char(M_.endo_names, 'y');
M_.endo_names_tex = char(M_.endo_names_tex, 'y');
M_.endo_names = char(M_.endo_names, 'i');
M_.endo_names_tex = char(M_.endo_names_tex, 'i');
M_.endo_names = char(M_.endo_names, 'k');
M_.endo_names_tex = char(M_.endo_names_tex, 'k');
M_.endo_names = char(M_.endo_names, 'a');
M_.endo_names_tex = char(M_.endo_names_tex, 'a');
M_.endo_names = char(M_.endo_names, 'lambda');
M_.endo_names_tex = char(M_.endo_names_tex, 'lambda');
M_.endo_names = char(M_.endo_names, 'tb_y');
M_.endo_names_tex = char(M_.endo_names_tex, 'tb\_y');
M_.endo_names = char(M_.endo_names, 'ca_y');
M_.endo_names_tex = char(M_.endo_names_tex, 'ca\_y');
M_.endo_names = char(M_.endo_names, 'r');
M_.endo_names_tex = char(M_.endo_names_tex, 'r');
M_.endo_names = char(M_.endo_names, 'util');
M_.endo_names_tex = char(M_.endo_names_tex, 'util');
M_.endo_names = char(M_.endo_names, 'beta_fun');
M_.endo_names_tex = char(M_.endo_names_tex, 'beta\_fun');
M_.endo_names = char(M_.endo_names, 'eta');
M_.endo_names_tex = char(M_.endo_names_tex, 'eta');
M_.param_names = 'gamma';
M_.param_names_tex = 'gamma';
M_.param_names = char(M_.param_names, 'omega');
M_.param_names_tex = char(M_.param_names_tex, 'omega');
M_.param_names = char(M_.param_names, 'rho');
M_.param_names_tex = char(M_.param_names_tex, 'rho');
M_.param_names = char(M_.param_names, 'sigmae');
M_.param_names_tex = char(M_.param_names_tex, 'sigmae');
M_.param_names = char(M_.param_names, 'delta');
M_.param_names_tex = char(M_.param_names_tex, 'delta');
M_.param_names = char(M_.param_names, 'psi_1');
M_.param_names_tex = char(M_.param_names_tex, 'psi\_1');
M_.param_names = char(M_.param_names, 'alpha');
M_.param_names_tex = char(M_.param_names_tex, 'alpha');
M_.param_names = char(M_.param_names, 'phi');
M_.param_names_tex = char(M_.param_names_tex, 'phi');
M_.param_names = char(M_.param_names, 'beta');
M_.param_names_tex = char(M_.param_names_tex, 'beta');
M_.param_names = char(M_.param_names, 'r_w');
M_.param_names_tex = char(M_.param_names_tex, 'r\_w');
M_.param_names = char(M_.param_names, 'd_bar');
M_.param_names_tex = char(M_.param_names_tex, 'd\_bar');
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.endo_nbr = 14;
M_.param_nbr = 11;
M_.orig_endo_nbr = 14;
M_.aux_vars = [];
M_.Sigma_e = zeros(1, 1);
M_.H = 0;
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
erase_compiled_function('idftest_static');
erase_compiled_function('idftest_dynamic');
M_.lead_lag_incidence = [
 1 5 0;
 0 6 0;
 0 7 0;
 0 8 19;
 0 9 0;
 2 10 20;
 3 11 0;
 0 12 21;
 0 13 0;
 0 14 0;
 4 15 0;
 0 16 22;
 0 17 23;
 0 18 24;]';
M_.nstatic = 5;
M_.nfwrd   = 5;
M_.npred   = 3;
M_.nboth   = 1;
M_.equations_tags = {
};
M_.exo_names_orig_ord = [1:1];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(14, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = NaN(11, 1);
M_.NNZDerivatives = zeros(3, 1);
M_.NNZDerivatives(1) = 59;
M_.NNZDerivatives(2) = 114;
M_.NNZDerivatives(3) = 283;
M_.params( 7 ) = 0.32;
alpha = M_.params( 7 );
M_.params( 3 ) = 0.86;
rho = M_.params( 3 );
M_.params( 8 ) = 0.0359;
phi = M_.params( 8 );
M_.params( 10 ) = 0.025;
r_w = M_.params( 10 );
M_.params( 1 ) = 2;
gamma = M_.params( 1 );
M_.params( 2 ) = 1.6;
omega = M_.params( 2 );
M_.params( 5 ) = 0.03;
delta = M_.params( 5 );
M_.params( 6 ) = 0.0635;
psi_1 = M_.params( 6 );
M_.params( 4 ) = 0.0201;
sigmae = M_.params( 4 );
M_.params( 9 ) = 1/(1+M_.params(10));
beta = M_.params( 9 );
h_ss   = ((1-alpha)*(alpha/(r_w+delta))^(alpha/(1-alpha)))^(1/(omega-1)); 
k_ss   = h_ss/(((r_w+delta)/alpha)^(1/(1-alpha)));
i_ss   = delta*k_ss;                                                     
y_ss   = (k_ss^alpha)*(h_ss^(1-alpha));                                   
c_ss = (1+r_w)^(1/ psi_1) + h_ss^ omega/ omega -1 ;
tb_ss  = y_ss-c_ss-i_ss;
M_.params( 11 ) = tb_ss/M_.params(10);
d_bar = M_.params( 11 );
d_ss   = d_bar;
%
% INITVAL instructions
%
options_.initval_file = 0;
oo_.steady_state( 11 ) = log((1-M_.params(9))/M_.params(9));
oo_.steady_state( 1 ) = d_ss;
oo_.steady_state( 3 ) = log(h_ss);
oo_.steady_state( 6 ) = log(k_ss);
oo_.steady_state( 4 ) = log(y_ss);
oo_.steady_state( 2 ) = log(c_ss);
oo_.steady_state( 5 ) = log(i_ss);
oo_.steady_state( 9 ) = 1-(exp(oo_.steady_state(2))+exp(oo_.steady_state(5)))/exp(oo_.steady_state(4));
oo_.steady_state( 12 ) = ((exp(oo_.steady_state(2))-M_.params(2)^(-1)*exp(oo_.steady_state(3))^M_.params(2))^(1-M_.params(1))-1)/(1-M_.params(1));
oo_.steady_state( 13 ) = (1+exp(oo_.steady_state(2))-M_.params(2)^(-1)*exp(oo_.steady_state(3))^M_.params(2))^(-M_.params(6));
oo_.steady_state( 14 ) = (-oo_.steady_state(12))/(1-oo_.steady_state(13));
oo_.steady_state( 8 ) = log((exp(oo_.steady_state(2))-exp(oo_.steady_state(3))^M_.params(2)/M_.params(2))^(-M_.params(1))-oo_.steady_state(14)*(-M_.params(6))*(1+exp(oo_.steady_state(2))-M_.params(2)^(-1)*exp(oo_.steady_state(3))^M_.params(2))^((-M_.params(6))-1));
if M_.exo_nbr > 0;
	oo_.exo_simul = [ones(M_.maximum_lag,1)*oo_.exo_steady_state'];
end;
if M_.exo_det_nbr > 0;
	oo_.exo_det_simul = [ones(M_.maximum_lag,1)*oo_.exo_det_steady_state'];
end;
oo_.dr.eigval = check(M_,options_,oo_);
steady;
%
% SHOCKS instructions
%
make_ex_;
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = (1)^2;
M_.sigma_e_is_diagonal = 1;
options_.drop = 1000;
options_.irf = 0;
options_.order = 3;
options_.periods = 50000;
var_list_=[];
info = stoch_simul(var_list_);
simulated_moments(M_,oo_, options_, var_list_, 3,  'fernandez-villaverde_et_al');
save('idftest_results.mat', 'oo_', 'M_', 'options_');


disp(['Total computing time : ' dynsec2hms(toc) ]);
diary off
