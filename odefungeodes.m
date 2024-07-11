function dydt = odefungeodes(t,y)
% t = le temps , vecteur


% Symboles de Christoffel calculés dans Christof.m
Gamma231=0.5; Gamma321=-0.5;
Gamma132=-0.5; Gamma312=0.5;
Gamma123=0.5; Gamma213=-0.5;
Gamma264=1; Gamma354=-1;
Gamma165=-1; Gamma345=1;
Gamma156=1; Gamma246=-1;

dydt = zeros(12,1);

% Equation to solve : dans se(3) (on exprime tout dans se(3))

% gamma^dotdot = - Gamma_ij^k gamma^dot i gamma^dot j
dydt(1)=-(Gamma231+Gamma321)*y(2)*y(3);
dydt(2)=-(Gamma132+Gamma312)*y(1)*y(3);
dydt(3)=-(Gamma123+Gamma213)*y(1)*y(2);
dydt(4)=-Gamma264*y(2)*y(6)-Gamma354*y(3)*y(5);
dydt(5)=-Gamma165*y(1)*y(6)-Gamma345*y(3)*y(4);
dydt(6)=-Gamma156*y(1)*y(5)-Gamma246*y(2)*y(4);

% gamma^dot = gamma^dot
dydt(7)=  y(1);
dydt(8)=  y(2);
dydt(9)=  y(3);
dydt(10)= y(4);
dydt(11)= y(5);
dydt(12)= y(6);

end