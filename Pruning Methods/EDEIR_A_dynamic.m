function [residual, g1, g2, g3] = EDEIR_A_dynamic(y, x, params, steady_state, it_)
%
% Status : Computes dynamic model for Dynare
%
% Inputs :
%   y         [#dynamic variables by 1] double    vector of endogenous variables in the order stored
%                                                 in M_.lead_lag_incidence; see the Manual
%   x         [nperiods by M_.exo_nbr] double     matrix of exogenous variables (in declaration order)
%                                                 for all simulation periods
%   steady_state  [M_.endo_nbr by 1] double       vector of steady state values
%   params    [M_.param_nbr by 1] double          vector of parameter values in declaration order
%   it_       scalar double                       time period for exogenous variables for which to evaluate the model
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the dynamic model equations in order of 
%                                          declaration of the equations.
%                                          Dynare may prepend auxiliary equations, see M_.aux_vars
%   g1        [M_.endo_nbr by #dynamic variables] double    Jacobian matrix of the dynamic model equations;
%                                                           rows: equations in order of declaration
%                                                           columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%   g2        [M_.endo_nbr by (#dynamic variables)^2] double   Hessian matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%   g3        [M_.endo_nbr by (#dynamic variables)^3] double   Third order derivative matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(12, 1);
T38 = exp(y(11))*exp(y(2))^params(7);
T42 = exp(y(7))^(1-params(7));
T62 = exp(y(7))^params(2);
T67 = (exp(y(6))-T62/params(2))^(-params(1));
T89 = params(9)*exp(y(19))*(1+params(7)*exp(y(17))/exp(y(10))-params(5)+params(8)*(exp(y(18))-exp(y(10))*params(5)));
lhs =y(5);
rhs =(1+exp(y(4)))*y(1)-exp(y(8))+exp(y(6))+exp(y(9))+params(8)/2*(exp(y(10))-exp(y(2)))^2;
residual(1)= lhs-rhs;
lhs =exp(y(8));
rhs =T38*T42;
residual(2)= lhs-rhs;
lhs =exp(y(10));
rhs =exp(y(9))+exp(y(2))*(1-params(5));
residual(3)= lhs-rhs;
lhs =exp(y(12));
rhs =params(9)*(1+exp(y(16)))*exp(y(19));
residual(4)= lhs-rhs;
lhs =T67;
rhs =exp(y(12));
residual(5)= lhs-rhs;
lhs =T62*T67;
rhs =exp(y(8))*(1-params(7))*exp(y(12));
residual(6)= lhs-rhs;
lhs =exp(y(12))*(1+params(8)*(exp(y(10))-exp(y(2))));
rhs =T89;
residual(7)= lhs-rhs;
lhs =y(11);
rhs =params(3)*y(3)+params(4)*x(it_, 1);
residual(8)= lhs-rhs;
lhs =y(13);
rhs =1-(exp(y(6))+exp(y(9)))/exp(y(8));
residual(9)= lhs-rhs;
lhs =y(14);
rhs =1/exp(y(8))*(y(5)-y(1));
residual(10)= lhs-rhs;
lhs =y(15);
rhs =params(6)*(exp(y(5)-params(11))-1);
residual(11)= lhs-rhs;
lhs =exp(y(16));
rhs =y(15)+params(10);
residual(12)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(12, 20);

  %
  % Jacobian matrix
  %

T125 = getPowerDeriv(exp(y(6))-T62/params(2),(-params(1)),1);
T135 = exp(y(7))*getPowerDeriv(exp(y(7)),params(2),1);
  g1(1,1)=(-(1+exp(y(4))));
  g1(1,5)=1;
  g1(1,6)=(-exp(y(6)));
  g1(1,8)=exp(y(8));
  g1(1,9)=(-exp(y(9)));
  g1(1,2)=(-(params(8)/2*(-exp(y(2)))*2*(exp(y(10))-exp(y(2)))));
  g1(1,10)=(-(params(8)/2*exp(y(10))*2*(exp(y(10))-exp(y(2)))));
  g1(1,4)=(-(exp(y(4))*y(1)));
  g1(2,7)=(-(T38*exp(y(7))*getPowerDeriv(exp(y(7)),1-params(7),1)));
  g1(2,8)=exp(y(8));
  g1(2,2)=(-(T42*exp(y(11))*exp(y(2))*getPowerDeriv(exp(y(2)),params(7),1)));
  g1(2,11)=(-(T38*T42));
  g1(3,9)=(-exp(y(9)));
  g1(3,2)=(-(exp(y(2))*(1-params(5))));
  g1(3,10)=exp(y(10));
  g1(4,12)=exp(y(12));
  g1(4,19)=(-(params(9)*(1+exp(y(16)))*exp(y(19))));
  g1(4,16)=(-(exp(y(19))*params(9)*exp(y(16))));
  g1(5,6)=exp(y(6))*T125;
  g1(5,7)=T125*(-(T135/params(2)));
  g1(5,12)=(-exp(y(12)));
  g1(6,6)=T62*exp(y(6))*T125;
  g1(6,7)=T67*T135+T62*T125*(-(T135/params(2)));
  g1(6,8)=(-(exp(y(8))*(1-params(7))*exp(y(12))));
  g1(6,12)=(-(exp(y(8))*(1-params(7))*exp(y(12))));
  g1(7,17)=(-(params(9)*exp(y(19))*params(7)*exp(y(17))/exp(y(10))));
  g1(7,18)=(-(params(9)*exp(y(19))*params(8)*exp(y(18))));
  g1(7,2)=exp(y(12))*params(8)*(-exp(y(2)));
  g1(7,10)=exp(y(12))*params(8)*exp(y(10))-params(9)*exp(y(19))*((-(exp(y(10))*params(7)*exp(y(17))))/(exp(y(10))*exp(y(10)))+params(8)*(-(exp(y(10))*params(5))));
  g1(7,12)=exp(y(12))*(1+params(8)*(exp(y(10))-exp(y(2))));
  g1(7,19)=(-T89);
  g1(8,3)=(-params(3));
  g1(8,11)=1;
  g1(8,20)=(-params(4));
  g1(9,6)=exp(y(6))/exp(y(8));
  g1(9,8)=(-(exp(y(8))*(exp(y(6))+exp(y(9)))))/(exp(y(8))*exp(y(8)));
  g1(9,9)=exp(y(9))/exp(y(8));
  g1(9,13)=1;
  g1(10,1)=1/exp(y(8));
  g1(10,5)=(-(1/exp(y(8))));
  g1(10,8)=(-((y(5)-y(1))*(-exp(y(8)))/(exp(y(8))*exp(y(8)))));
  g1(10,14)=1;
  g1(11,5)=(-(params(6)*exp(y(5)-params(11))));
  g1(11,15)=1;
  g1(12,15)=(-1);
  g1(12,16)=exp(y(16));

if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],12,400);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],12,8000);
end
end
end
end
