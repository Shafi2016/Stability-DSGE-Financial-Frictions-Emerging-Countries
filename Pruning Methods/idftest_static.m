function [residual, g1, g2] = idftest_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 14, 1);

%
% Model equations
%

T36 = exp(y(7))*exp(y(6))^params(7)*exp(y(3))^(1-params(7));
T52 = exp(y(2))-exp(y(3))^params(2)/params(2);
T55 = T52^(-params(1));
T61 = exp(y(3))^params(2)*params(2)^(-1);
T62 = 1+exp(y(2))-T61;
T65 = (-params(6))*T62^((-params(6))-1);
T75 = exp(y(3))^(params(2)-1);
T90 = exp(y(8))*y(13)*(1+exp(y(4))*params(7)/exp(y(6))-params(5));
T148 = getPowerDeriv(T52,(-params(1)),1)*(-(exp(y(3))*getPowerDeriv(exp(y(3)),params(2),1)/params(2)));
T150 = (-(params(2)^(-1)*exp(y(3))*getPowerDeriv(exp(y(3)),params(2),1)));
T152 = (-params(6))*getPowerDeriv(T62,(-params(6))-1,1)*T150;
T178 = (-(exp(y(4))*(1-params(7))*exp(y(8))/exp(y(3))));
lhs =y(1);
rhs =y(1)*(1+exp(y(11)))-exp(y(4))+exp(y(2))+exp(y(5));
residual(1)= lhs-rhs;
lhs =exp(y(4));
rhs =T36;
residual(2)= lhs-rhs;
lhs =exp(y(6));
rhs =exp(y(5))+exp(y(6))*(1-params(5));
residual(3)= lhs-rhs;
lhs =exp(y(8));
rhs =exp(y(8))*(1+exp(y(11)))*y(13);
residual(4)= lhs-rhs;
lhs =exp(y(8));
rhs =T55-y(14)*T65;
residual(5)= lhs-rhs;
lhs =y(14);
rhs =(-y(12))+y(13)*y(14);
residual(6)= lhs-rhs;
lhs =T55*T75+y(14)*T65*(-T75);
rhs =exp(y(4))*(1-params(7))*exp(y(8))/exp(y(3));
residual(7)= lhs-rhs;
lhs =exp(y(8));
rhs =T90;
residual(8)= lhs-rhs;
lhs =y(7);
rhs =y(7)*params(3)+params(4)*x(1);
residual(9)= lhs-rhs;
lhs =y(13);
rhs =T62^(-params(6));
residual(10)= lhs-rhs;
lhs =y(12);
rhs =((exp(y(2))-T61)^(1-params(1))-1)/(1-params(1));
residual(11)= lhs-rhs;
lhs =exp(y(11));
rhs =params(10);
residual(12)= lhs-rhs;
lhs =y(9);
rhs =1-(exp(y(2))+exp(y(5)))/exp(y(4));
residual(13)= lhs-rhs;
residual(14) = y(10);
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(14, 14);

  %
  % Jacobian matrix
  %

  g1(1,1)=1-(1+exp(y(11)));
  g1(1,2)=(-exp(y(2)));
  g1(1,4)=exp(y(4));
  g1(1,5)=(-exp(y(5)));
  g1(1,11)=(-(y(1)*exp(y(11))));
  g1(2,3)=(-(exp(y(7))*exp(y(6))^params(7)*exp(y(3))*getPowerDeriv(exp(y(3)),1-params(7),1)));
  g1(2,4)=exp(y(4));
  g1(2,6)=(-(exp(y(3))^(1-params(7))*exp(y(7))*exp(y(6))*getPowerDeriv(exp(y(6)),params(7),1)));
  g1(2,7)=(-T36);
  g1(3,5)=(-exp(y(5)));
  g1(3,6)=exp(y(6))-exp(y(6))*(1-params(5));
  g1(4,8)=exp(y(8))-exp(y(8))*(1+exp(y(11)))*y(13);
  g1(4,11)=(-(exp(y(8))*exp(y(11))*y(13)));
  g1(4,13)=(-((1+exp(y(11)))*exp(y(8))));
  g1(5,2)=(-(exp(y(2))*getPowerDeriv(T52,(-params(1)),1)-y(14)*(-params(6))*exp(y(2))*getPowerDeriv(T62,(-params(6))-1,1)));
  g1(5,3)=(-(T148-y(14)*T152));
  g1(5,8)=exp(y(8));
  g1(5,14)=T65;
  g1(6,12)=1;
  g1(6,13)=(-y(14));
  g1(6,14)=1-y(13);
  g1(7,2)=T75*exp(y(2))*getPowerDeriv(T52,(-params(1)),1)+y(14)*(-T75)*(-params(6))*exp(y(2))*getPowerDeriv(T62,(-params(6))-1,1);
  g1(7,3)=T75*T148+T55*exp(y(3))*getPowerDeriv(exp(y(3)),params(2)-1,1)+y(14)*((-T75)*T152+T65*(-(exp(y(3))*getPowerDeriv(exp(y(3)),params(2)-1,1))))-(-(exp(y(3))*exp(y(4))*(1-params(7))*exp(y(8))))/(exp(y(3))*exp(y(3)));
  g1(7,4)=T178;
  g1(7,8)=T178;
  g1(7,14)=T65*(-T75);
  g1(8,4)=(-(exp(y(8))*y(13)*exp(y(4))*params(7)/exp(y(6))));
  g1(8,6)=(-(exp(y(8))*y(13)*(-(exp(y(6))*exp(y(4))*params(7)))/(exp(y(6))*exp(y(6)))));
  g1(8,8)=exp(y(8))-T90;
  g1(8,13)=(-(exp(y(8))*(1+exp(y(4))*params(7)/exp(y(6))-params(5))));
  g1(9,7)=1-params(3);
  g1(10,2)=(-(exp(y(2))*getPowerDeriv(T62,(-params(6)),1)));
  g1(10,3)=(-(getPowerDeriv(T62,(-params(6)),1)*T150));
  g1(10,13)=1;
  g1(11,2)=(-(exp(y(2))*getPowerDeriv(exp(y(2))-T61,1-params(1),1)/(1-params(1))));
  g1(11,3)=(-(getPowerDeriv(exp(y(2))-T61,1-params(1),1)*T150/(1-params(1))));
  g1(11,12)=1;
  g1(12,11)=exp(y(11));
  g1(13,2)=exp(y(2))/exp(y(4));
  g1(13,4)=(-(exp(y(4))*(exp(y(2))+exp(y(5)))))/(exp(y(4))*exp(y(4)));
  g1(13,5)=exp(y(5))/exp(y(4));
  g1(13,9)=1;
  g1(14,10)=1;
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],14,196);
end
end
