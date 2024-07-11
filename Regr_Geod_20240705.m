% 02/07/2024 
% Régression géodésique sur SE(3) 
% MAIS on prend un point de départ x0 et les traj réelles et bruitées
% générées à partir de x0 (x_n = g_n * x_0) et non pas de proche en proche 
% 
clear all
close all
%% Manoeuvres et vitesses 
N_mes = 50; % Nombre de mesures (nombre de déplacements dans l'intervalle [0;1]) 
n_maneuv = 1; 

x0 = [1 ; 1 ; 1; 1]; % Point de départ du mobile 

% Déplacement à vitesse constante ( vitesse ramenée à l'identité constante )
omegax = 0.0; 
omegay = 2.0; 
omegaz = 0; 
vx = 0 ; 
vy = 1 ; 
vz = 0 ; 
% on génère le vecteur tangent dans se(3) correspondant au mouvement
ksi1 = hat_map_se3([omegax  omegay  omegaz vx vy vz]) ; 

% % Déplacement à vitesse constante ( vitesse ramenée à l'identité constante )
% omegax = 0;
% omegay = 0;
% omegaz = 0;
% vx = 20 ;
% vy = 0 ;
% vz = 0 ;
% % on génère le vecteur tangent dans se(3) correspondant au mouvement
% ksi2 = hat_map_se3([omegax  omegay  omegaz vx vy vz]) ;

%% Génération des données 

% Ecart-type du bruit additif selon une loi N(0, sigma^2) dans se(3) : 
sigma = 0.001; 

% Initialisation des données réelles et des données mesurées dans R^3 (le 1 complète la 4ème
% ligne) et des élements de SE(3) réels pour chaque observation
xreel = zeros(4,n_maneuv*N_mes); % Tableau 4-par-N_mes contenant les positions observées bruitées
greel = zeros(4,4,n_maneuv*N_mes);  % 1er élement d'un tableau de dimension 4 par 4 par N_mes
xreel(4,:) = 1; % Dernière ligne uniquement des 1
xreel(:,1) = x0; % 1ère colonne : point de départ
xmes = zeros(4,n_maneuv*N_mes); % Tableau 4-par-N_mes contenant les positions observées bruitées
gmes = zeros(4,4,n_maneuv*N_mes);  % 1er élement d'un tableau de dimension 4 par 4 par N_mes
xmes(4,:) = 1; % Dernière ligne uniquement des 1
xmes(:,1) = x0; % 1ère colonne : point de départ

% Génération des données réelles 
for i=1:N_mes
    greel(:,:,i) = expm(i/N_mes.*ksi1); % On crée l'élement du groupe SE(3) qui représente le déplacement entre les instants i et i+1 
    xreel(:,i+1) = greel(:,:,i) * x0; % Déplacement par action du groupe SE(3) sur R^3 x {1}
end
% for i=N_mes+1:2*N_mes
%     greel(:,:,i) = expm(1/N_mes.*ksi2); % On crée l'élement du groupe SE(3) qui représente le déplacement entre les instants i et i+1
%     xreel(:,i+1) = greel(:,:,i) * xreel(:,i); % Déplacement par action du groupe SE(3) sur R^3 x {1}
% end


% Génération des données bruitées 
for i=1:N_mes 
    b = sigma.* [hat_map_so3(randn(1,3)) , randn(3,1); 0 0 0 0]; % Génération du bruit b dans se(3) sur les 6 composantes
    gmes(:,:,i) = expm(i/N_mes.*ksi1 + b); % On crée l'élement du groupe SE(3) qui représente le déplacement entre les instants i et i+1 
    xmes(:,i+1) = gmes(:,:,i) * x0; % Déplacement par action du groupe SE(3) sur R^3 x {1}
end
% for i=N_mes+1:2*N_mes
%     b = sigma.*[hat_map(randn(1,3)), [randn(3,1)] ; 0 0 0 0]; % Génération du bruit b dans se(3) sur les 6 composantes
%     gmes(:,:,i) = expm(1/N_mes.*ksi2 + b); % On crée l'élement du groupe SE(3) qui représente le déplacement entre les instants i et i+1
%     xmes(:,i+1) = gmes(:,:,i) * xreel(:,i); % Déplacement par action du groupe SE(3) sur R^3 x {1}
% end
% figure
% plot(xreel(1,:),xreel( 2,:),'+-');
% hold on;
% plot(xmes(1,:), xmes(2,:),'r+');
% figure
% plot3(xreel(1,:),xreel( 2,:),xreel(3,:),'+-');
% hold on;
% plot3(xmes(1,:), xmes(2,:),xmes(3,:),'r+'); 
% clc 

%% Algorithme pour déterminer la géodésique 
n=[]; 
tic 
%for k = [1] 
    
N_calc = 10; % pas élémentaire pour la résolution de l'ED 
epsilon = 1e-3; % pas de descente de gradient 
Nit = 1300; % nombre d'itérations 

% Initialisations 

n=[]; 
M=0; 
n0=[]; 
n1=[]; 
r=[]; 
ner=[]; 
xest_R3 = zeros(4,n_maneuv*N_mes+1); 
xest_R3(4,:) = 1; 
xest_R3(:,1) = x0; 

% Calcul de la moyenne de Fréchet des points SE(3) observés (mesurés) 
Nsteps=5; 
frechet_mean = eye(4); 
    for k=1:Nsteps 
        delta_i = zeros(4,4,N_mes); 
            for i=1:N_mes  
                delta_i(:,:,i) = logm(inv(frechet_mean)* gmes(:,:,i)); 
            end
        d = cumsum( delta_i ,3 ); 
        delta = 1/N_mes.*d(:,:,N_mes); 
        frechet_mean = frechet_mean*expm(delta); 
    end 

gamma_est_SE3 = zeros(4,4,n_maneuv*N_mes); 

gamma_dot_0 = zeros(4,4); 

%gamma_0 = eye(4,4); 
gamma_0 = frechet_mean; 
%gamma_0 = gmes(:,:,floor(N_mes*rand));


for l=1:Nit 
    if(mod(l,10)== 0)
        l
    end
    % Etape 1 : calcul de la géodésique 
    % Résolution équations des géodésiques pour obtenir gamma_e_dot(t) (\in algèbre de Lie) et 
    % gamma_e(t) (\in algèbre de Lie) pour tout t sur [0;1] (contrainte gamma(t) = géodésique) 
    
    % Initlisation de gamma_e_dot(0) et gamma_e(0) 
    if (l==1)
            gamma_dot_e_0_vect = vee_map_se3(inv(gamma_0)*gamma_dot_0); % 1) on ramène gamma_dot_0 dans se(3) 
            gamma_dot_e_0_vect = gamma_dot_e_0_vect';
    end
   
    %gamma_dot_0 = gamma_0 * hat_map_se3(gamma_dot_e_0_vect)  % 3) on renvoie dans T_gamma(0)_SE(3) (espace tgt à gamma(0)) 
    %gamma_e_0 = SE3_se3_back(gamma_0)'; 
    
    [ttt,yyy] = ode45(@(t,y) odefungeodes(t,y) , linspace(0,1,n_maneuv*N_mes+1), [gamma_dot_e_0_vect' ; zeros(6,1) ]); 
    % yyy(1 à 6)  gamma_e_dot 
    % yyy(7 à 12) gamma_e 
    gamma_e = yyy(:,7:12); 
    gamma_dot_e = yyy(:,1:6); 
 
    for i=1:n_maneuv*N_mes+1 
        gamma_est_SE3(:,:,i) = gamma_0 * expm(hat_map_se3(gamma_e(i,:))); % gamma = expm (gamma_e) avec gamma_e dans se(3) 
        xest_R3(:,i) = gamma_est_SE3(:,:,i) * x0; 
    end 
    
    %l = 1 Ok vérifié 
    
    % Etape 2 : transport de la forme de t=1 vers t=0 (sens inverse) 
    % Minimisation du critère de régression 
    tt = []; 
    yy = []; 
    % y0 condition initale de l'EDO  application de lambda0(1-) = 0 + log
    % (gamma(1) , mesure à l'instant 1) ramené dans se3 
    temp = inv(gamma_est_SE3(:,:,n_maneuv*N_mes+1))*riemlog_se3(gamma_est_SE3(:,:,n_maneuv*N_mes+1) , gmes(:,:,n_maneuv*N_mes));
    y0 = [ vee_map_se3(temp)' , zeros(1,6)]'; 
    %y0 = [ zeros(1,6) , zeros(1,6)]'; 
    
    for j=1:N_mes 
   
        tspan = linspace( 1 - (j-1)/N_mes , 1 - j/N_mes , N_calc); % on résout sur [t (j+1); t(j)] (temps en sens inverse)
        
        [t,y] = ode45(@(t,y) odefunlambda(t,y, gamma_dot_e(N_mes-j+1,:) ) , tspan, y0);
        %[t,y] = ode45(@(t,y) odefunlambda(t,y, gamma_dot_e(N_mes-j+1,:) ) , [1 - (j-1)/N_mes ; 1 - j/N_mes], y0);
        tt = [tt ; 1 - (j-1)/N_mes];
        yy = [yy ; y(end,:)] ;
        if (j ~= N_mes) 
            % y0 condition initale de l'EDO application de lambda0(tn -) = lambda0(tn +) + [log (gamma(tn) , mesure à l'instant tn)] ramené dans TeG = g )
            temp_e = vee_map_se3( inv(gamma_est_SE3(:,:,n_maneuv*N_mes-j)) * riemlog_se3(gamma_est_SE3(:,:,n_maneuv*N_mes-j) , gmes(:,:,n_maneuv*N_mes-j)) )'; 
            y0 = yy(end,:)' +  [ temp_e , zeros(1,6)]' ;
            %y0 = y(end,:)' +  [ temp_e , zeros(1,6)]' ;
        end
        
    end 
    % On a transporté la forme en t=0 
    
    % Mise à jour de gamma(0) et gamma_dot(0) 
    lambda0_0_vect = yy(end,1:6); % variable adjointe en t=0 qui sert à mettre à jour gamma_0 (dans SE(3))
    lambda0_0_mat = hat_map_se3(lambda0_0_vect); 
    lambda1_0_vect = yy(end,7:12); % variable adjointe en t=0 qui sert à mettre à jour X = gamma'_e(0) (dans se(3))
    lambda1_0_mat = hat_map_se3(lambda1_0_vect); 
    gamma_0 = gamma_0 * expm(epsilon.*lambda0_0_mat); % mise à jour de gamma_0 (dans SE(3) ) 
    gamma_dot_e_0_vect = gamma_dot_e(1,:); 
    gamma_dot_e_0_vect = gamma_dot_e_0_vect + epsilon .* lambda1_0_vect; % 2) on met à jour gamma_dot_0 dans se(3) 
    gamma_dot_0 = gamma_0 * hat_map_se3(gamma_dot_e_0_vect);  % 3) on renvoie dans T_gamma(0)_SE(3) (espace tgt à gamma(0))
    
    n0 = [n0 ; norm(lambda0_0_vect)]; 
    n1 = [n1 ; norm(lambda1_0_vect)]; 
    
    if (mod(l,100) == 0) % toutes les 10 itérations ... 
        epsilon=epsilon; % ... on diminue le pas de descente (meilleure convergence)
    end
    
    r=[];
    for i=1:n_maneuv*N_mes 
        r = [r ; norm( riemlog_se3( gamma_est_SE3(:,:,i+1) , gmes(:,:,i) ) ) ]; % écart au carré entre gamma(tj) et mesure j
    end 
    ner = [ner;sum(r)]; % norme de l'erreur Riemannienne 

end

%% Après la résolution 
%clc 
figure 
plot([1:Nit],log(ner)); 
title('log de Riemannian error norm') 
xlabel('Iterations') 
% figure
% plot([1:Nit],n1/N_mes); 
% xlabel('Iterations') 


% xest = zeros(4,n_maneuv*N_mes); 
% xest(4,:) = 1; 
% xest(:,1) = x0; 
% 
% % Résolution de l'équation des géodésiques cf fichier odefungeodes  pour
% % obtenir gamm_e_dot et gamma_e 
% % yyy(1 à 6) = gamma_e_dot 
% % yyy(7 à 12)= gamma_e 
% [ttt,yyy] = ode45(@(t,y) odefungeodes(t,y) , linspace(0,1,n_maneuv*N_mes), [gamma_dot_0' ; SE3_se3_back(gamma_0)']); 
% 
% gamma_e = yyy(:,7:12); 
% gamma_mat = zeros(4,4,n_maneuv*N_mes); 
% for i=1:n_maneuv*N_mes 
%     gamma_mat(:,:,i) = expm(1/N_mes.*hat_map_se3(gamma_e(i,:))); % gamma = expm (gamma_e) avec gamma_e dans se(3) 
%     % gamma_mat(:,:,i) = expm(hat_map_se3(gamma_0)); % gamma = expm (gamma_e) avec gamma_e dans se(3) 
%     xest(:,i+1) = gamma_mat(:,:,i) * xest(:,i); 
% end 

figure 
plot(xest_R3(1,:),xest_R3(2,:),'g-o'); 
hold on; 
plot(xmes(1,2:end), xmes(2,2:end),'r+'); 
hold on; 
plot(xreel(1,2:end), xreel(2,2:end),'b+-'); 
legend('Estimate','Measured','Actual') 

figure 
plot3(xest_R3(1,:),xest_R3( 2,:),xest_R3(3,:),'g-o'); 
hold on; 
plot3(xmes(1,:), xmes(2,:),xmes(3,:),'r+'); 
hold on; 
plot3(xreel(1,:), xreel(2,:),xreel(3,:),'b+-' ); 
%axis([-30 30 -30 30 -10 10]) 
legend('Estimate','Measured','Actual') 
n=[n, (norm(xmes'-mean(xmes'))-norm(xmes'-xest_R3')) / norm(xmes'-mean(xmes'))]; 
[M, indice] = max(n); 
%M 

%% Calcul du R^2 Riemannien et Euclidien 

frechet_var = 0; 
for i=1:N_mes 
    delta_i(:,:,i) = logm( (inv(frechet_mean)) * gmes(:,:,i) ); 
    frechet_var = frechet_var + trace(delta_i(:,:,i)'*delta_i(:,:,i)); 
end 
frechet_var = frechet_var / N_mes; 

num = 0; 
for i=1:N_mes 
    delta_i(:,:,i) = riemlog_se3( gamma_est_SE3(:,:,i+1), gmes(:,:,i) ); 
    num = num + trace(delta_i(:,:,i)'*delta_i(:,:,i)); 
end 
num = num / N_mes; 

%mu_frechet_hat 

r_squared = 1 - num / frechet_var 

toc  



%%
% writematrix(xest,'yest_R3.csv','Delimiter','space'); 
% writematrix(xmes,'ymes_R3.csv','Delimiter','space'); 
% writematrix(gamma_mat,'yest_SE3.csv','Delimiter','space'); 
% writematrix(gmes,'ymes_SE3.csv','Delimiter','space'); 