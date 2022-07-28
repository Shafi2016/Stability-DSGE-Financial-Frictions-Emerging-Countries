function [residual, g1, g2, g3] = EDEIR_A_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Inputs : 
%   y         [M_.endo_nbr by 1] double    vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1] double     vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1] double   vector of parameter values in declaration order
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the static model equations 
%                                          in order of declaration of the equations.
%                                          Dynare may prepend or append auxiliary equations, see M_.aux_vars
%   g1        [M_.endo_nbr by M_.endo_nbr] double    Jacobian matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%   g2        [M_.endo_nbr by (M_.endo_nbr)^2] double   Hessian matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%   g3        [M_.endo_nbr by (M_.endo_nbr)^3] double   Third derivatives matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 12, 1);

%
% Model equations
%

T31 = exp(y(7))*exp(y(6))^params(7);
T35 = exp(y(3))^(1-params(7));
T50 = exp(y(3))^params(2);
T55 = (exp(y(2))-T50/params(2))^(-params(1));
T70 = exp(y(8))*params(9)*(1+exp(y(4))*params(7)/exp(y(6))-params(5)+params(8)*(exp(y(5))-exp(y(6))*params(5)));
lhs =y(1);
rhs =y(1)*(1+exp(y(12)))-exp(y(4))+exp(y(2))+exp(y(5));
residual(1)= lhs-rhs;
lhs =exp(y(4));
rhs =T31*T35;
residual(2)= lhs-rhs;
lhs =exp(y(6));
rhs =exp(y(5))+exp(y(6))*(1-params(5));
residual(3)= lhs-rhs;
lhs =exp(y(8));
rhs =exp(y(8))*(1+exp(y(12)))*params(9);
residual(4)= lhs-rhs;
lhs =T55;
rhs =exp(y(8));
residual(5)= lhs-rhs;
lhs =T50*T55;
rhs =exp(y(4))*(1-params(7))*exp(y(8));
residual(6)= lhs-rhs;
lhs =exp(y(8));
rhs =T70;
residual(7)= lhs-rhs;
lhs =y(7);
rhs =y(7)*params(3)+params(4)*x(1);
residual(8)= lhs-rhs;
lhs =y(9);
rhs =1-(exp(y(2))+exp(y(5)))/exp(y(4));
residual(9)= lhs-rhs;
residual(10) = y(10);
lhs =y(11);
rhs =params(6)*(exp(y(1)-params(11))-1);
residual(11)= lhs-rhs;
lhs =exp(y(12));
rhs =y(11)+params(10);
residual(12)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(12, 12);

  %
  % Jacobian matrix
  %

T102 = getPowerDeriv(exp(y(2))-T50/params(2),(-params(1)),1);
T112 = exp(y(3))*getPowerDeriv(exp(y(3)),params(2),1);
  g1(1,1)=1-(1+exp(y(12)));
  g1(1,2)=(-exp(y(2)));
  g1(1,4)=exp(y(4));
  g1(1,5)=(-exp(y(5)));
  g1(1,12)=(-(y(1)*exp(y(12))));
  g1(2,3)=(-(T31*exp(y(3))*getPowerDeriv(exp(y(3)),1-params(7),1)));
  g1(2,4)=exp(y(4));
  g1(2,6)=(-(T35*exp(y(7))*exp(y(6))*getPowerDeriv(exp(y(6)),params(7),1)));
  g1(2,7)=(-(T31*T35));
  g1(3,5)=(-exp(y(5)));
  g1(3,6)=exp(y(6))-exp(y(6))*(1-params(5));
  g1(4,8)=exp(y(8))-exp(y(8))*(1+exp(y(12)))*params(9);
  g1(4,12)=(-(exp(y(8))*exp(y(12))*params(9)));
  g1(5,2)=exp(y(2))*T102;
  g1(5,3)=T102*(-(T112/params(2)));
  g1(5,8)=(-exp(y(8)));
  g1(6,2)=T50*exp(y(2))*T102;
  g1(6,3)=T55*T112+T50*T102*(-(T112/params(2)));
  g1(6,4)=(-(exp(y(4))*(1-params(7))*exp(y(8))));
  g1(6,8)=(-(exp(y(4))*(1-params(7))*exp(y(8))));
  g1(7,4)=(-(exp(y(8))*params(9)*exp(y(4))*params(7)/exp(y(6))));
  g1(7,5)=(-(exp(y(8))*params(9)*exp(y(5))*params(8)));
  g1(7,6)=(-(exp(y(8))*params(9)*((-(exp(y(6))*exp(y(4))*params(7)))/(exp(y(6))*exp(y(6)))+params(8)*(-(exp(y(6))*params(5))))));
  g1(7,8)=exp(y(8))-T70;
  g1(8,7)=1-params(3);
  g1(9,2)=exp(y(2))/exp(y(4));
  g1(9,4)=(-(exp(y(4))*(exp(y(2))+exp(y(5)))))/(exp(y(4))*exp(y(4)));
  g1(9,5)=exp(y(5))/exp(y(4));
  g1(9,9)=1;
  g1(10,10)=1;
  g1(11,1)=(-(params(6)*exp(y(1)-params(11))));
  g1(11,11)=1;
  g1(12,11)=(-1);
  g1(12,12)=exp(y(12));
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],12,144);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],12,1728);
end
end
end
end
