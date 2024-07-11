function dydt = odefunlambda(t,y,X) 
% t = le temps , vecteur 
% y = le vecteur colonne de longueur 12 des composantes 
% X vecteur tangent (vitesse exprimée dans l'algèbre de Lie) 


% Vecteur y' à 12 composantes
dydt = zeros(12,1);

% Equation to solve : dans se(3) (on exprime tout dans se(3))

% Vecteur à droite
% y(1) à y(6) désignent lambda_0 1 à 6 
% y(7) à y(12) désignent lambda_1 1 à 6 

% Vecteur à gauche 
% dydt(1 à 6) désignent les lambda_0^dot 1 à 6 
% dydt(7) à dydt(12) : lambda_1^dot 1 à 6 
ad_X = [ hat_map_so3(X(1:3)) , zeros(3,3)  ; hat_map_so3(X(4:6)) , hat_map_so3(X(1:3)) ];
dydt(1:6) = - ad_X' * y(1:6) ;
ad_Y = [ hat_map_so3(y(1:3)) , zeros(3,3)  ; hat_map_so3(y(4:6)) , hat_map_so3(y(1:3)) ];
% lambda_1^dot = - lambda_0 + sym*_X lambda_1
%sym_etoile_X = ad_X' - ad_Y;
%sym_etoile_X = ad_X' * y(7:12) - ad_Y * X';
sym_etoile_X = -ad_X * y(7:12) + ad_Y' * X';
dydt(7:12) = - y(1:6) - sym_etoile_X ;
 
end 