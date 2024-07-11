function [R2, dist] = rSquaredSE3(data, fitted, meanTransform) % Calcul du R^2 Riemannien dans SE(3) 

frechet_mean = meanTransform; 

N_mes = size(data,3); 

delta_i = zeros(4,4,N_mes); 
var_tot = 0; 
for i=1:N_mes  
    %delta_i(:,:,i) = logm( (inv(data(:,:,i))) * frechet_mean  ); 
    %var_tot = var_tot + trace( delta_i(:,:,i)' * delta_i(:,:,i) ); 
    var_tot = var_tot + se3Distance( data(:,:,i) , frechet_mean )^2; 
end 
var_tot = var_tot / N_mes  % Calcul OK vérifié 28/06/24 


%delta_i = zeros(4,4,N_mes); 
var_expl = 0; 
for i=1:N_mes 
    % delta_i(:,:,i) = logm( inv(data(:,:,i)) * fitted(:,:,i) ); 
    % var_expl = var_expl + trace( delta_i(:,:,i)' * delta_i(:,:,i) ); 
    var_expl = var_expl + se3Distance( data(:,:,i) , fitted(:,:,i) )^2; 
end 
var_expl = var_expl / N_mes 

r_squared = 1 - var_expl / var_tot; 
R2 = r_squared; 
end 