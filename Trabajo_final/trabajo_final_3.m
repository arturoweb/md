%
% Trabajo final - Macroeconomía Dinámica
% Profesor: Carlos Urrutia
% ITAM, 2015
%
% Equipo:
% Omar Trejo, 119711
% Alejandro Cerecero, 000000
% Arturo Reynoso, 000000
%

%% Limpiar memoria

clear all
close all
clc

%% Parámetros del modelo

alpha      = 0.40;
beta       = 0.987;
delta      = 0.012;
gamma      = 0.64;
rho        = 0.95;
sigma_eps  = 0.0024;
sigma_zeta = 0.007;

%% Parámetros de la solución

max_iter   = 1000;
grid       = 500;
tol        = 1e-3;
T          = 100;
A          = 10;
kss        = ((A*beta*alpha) / ...
              (1 - (1 - delta)*beta))^(1/(1 - alpha));
k0         = (2/3)*kss;
p          = 500;  % Tamaño de la malla

%%%%%%%%%%%%%%%    Inciso (a)    %%%%%%%%%%%%%%%

% Cinco estados para el shock
theta = [0.0231 -0.0115 0.0115 0.0231];
q     = length(theta);

% Matriz de transición
Pi = [0.9727 0.0273      0      0      0;
      0.0041 0.9806 0.0153      0      0;
           0 0.0082 0.9837 0.0082      0;
           0      0 0.0153 0.9806 0.0041;
           0      0      0 0.0273 0.9727];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha    = 0.35;
beta     = 0.99;
delta    = 0.06;
A        = 10;
max_iter = 1000;
p        = 500;
tol      = 1e-2;
T        = 50;

theta = [0 0.5 1];
q     = length(theta);
Pi    = [0.5 0.3 0.2;
         0.1 0.7 0.2;
         0.1 0.4 0.5];

% Estado estacionario sin incertidumbre
kss  = ((A*beta*alpha)/(1 - (1 - delta)*beta))^(1/(1 - alpha));
k0   = kss;

% TODO: (otn) ¿Por qué este tamaño y no otro?
k    = linspace(0, 1.2*exp(1)*kss, p);

k(1) = 0.000001;

value = zeros(p*q, p, 3);
for i=1:p
    value((i - 1)*q + 1, :, 1) = k(i)*ones(1, p);
    value((i - 1)*q + 2, :, 1) = k(i)*ones(1, p);
    value((i - 1)*q + 3, :, 1) = k(i)*ones(1, p);
end

% Matriz de capitales por filas
value(:, :, 2) = ones(p*q, 1)*k;

vector = theta';
for i = 1:p-1
    vector = [vector; theta'];
end

% TODO: (otn) matriz de ...
value(:, :, 3) = vector*ones(1, p);

M      = zeros(p*q,p);
Mzeros = zeros(p*q,p);

% TODO: (otn) ¿qué es esto?
%     F(Sij,Kl)=F(Ki,Zj,Kl)=u(exp(Zj)*f(Ki)+(1-delta)Ki-Kl);

M = log(max(exp(value(:, :, 3)).*A.*value(:, :,1).^(alpha) - ...
            value(:, :, 2) + (1 - delta).*value(:, :, 1), ...
            0));

% Repetición de la identidad
I = eye(q);
E = I;
for i = 1:p-1
    E = [E; I];
end

% Algoritmo de optimización
V0         = zeros(p*q, 1);
finalizado = 0;
iter       = 1;

while finalizado == 0 && iter < max_iter
    [V1, G] = max((M + beta*(E*(Pi*reshape(V0, q, p))))');
    V1 = V1';
    G  = G';
    if norm(V0 - V1) < tol
        finalizado = 1;
    end
    disp(['Iter = ', num2str(iter), ' Error= ', num2str(norm(V0 - V1))])
    V0 = V1;
    iter = iter + 1;
end

G1 = G;
G  = zeros(p, q);
V  = zeros(p, q);

for i=1:p
    G(i,1) = G1((i - 1)*q + 1);
    G(i,2) = G1((i - 1)*q + 2);
    G(i,3) = G1((i - 1)*q + 3);
    V(i,1) = V1((i - 1)*q + 1);
    V(i,2) = V1((i - 1)*q + 2);
    V(i,3) = V1((i - 1)*q + 3);
end

G(1,1) = 1;
G(1,2) = 1;
G(1,3) = 1;
V(1,1) = V(2,1);
V(1,2) = V(2,2);
V(1,3) = V(2,3);

for i = 1:p
    c_1(i) = exp(theta(1))*A*k(i)^alpha - ...
             k(G(i, 1)) + (1 - delta)*k(i);
    c_2(i) = exp(theta(2))*A*k(i)^alpha - ...
             k(G(i, 2)) + (1 - delta)*k(i);
    c_3(i) = exp(theta(3))*A*k(i)^alpha - ...
             k(G(i, 3)) + (1 - delta)*k(i);
end

figure(1);

subplot(2, 1, 1);
plot(k, k(G(:, 1)), 'b', ...
     k, k(G(:, 2)), 'r', ...
     k, k(G(:, 3)), 'g');
     legend('theta = 0', 'theta = 0.5', 'theta = 1');
     title('Capital optimo');

subplot(2, 1, 2);
plot(k, c_1, 'b', ...
     k, c_2, 'r', ...
     k, c_3, 'g');
     legend('theta = 0', 'theta = 0.5', 'theta = 1');
     title('Consumo optimo');

%%%%%%%%%%%%%%%    Inciso (b)    %%%%%%%%%%%%%%%

% Nota: Ajustar, viene del inciso 5 del laboratorio 2.

% Distribución invariante
A = Pi^10000;

T_stat = 500;
Nsim   = 500;

theta = [1 2 3];
Pi0   = [0 1 0];
Pi    = [0.5 0.3 0.2;
         0.1 0.7 0.2;
         0.1 0.4 0.5];

for m = 1:Nsim
    [X_path, states] = markov(estado_choque, Pi, Pi0, T_stat);
    shocks_stat = X_path';

    kt_stat           = zeros(T_stat + 1, 1);
    ind_stat          = zeros(T_stat, 1);
    [aux,ind_stat(1)] = min(k < k0);
    kt_stat(1)        = k(ind_stat(1));

    for t = 1:T_stat
        ind_stat(t + 1) = G(ind_stat(t), shocks_stat(t));
        kt_stat(t + 1)  = k(ind_stat(t + 1));
    end

    it_stat = zeros(T_stat, 1);
    ct_stat = zeros(T_stat, 1);
    yt_stat = zeros(T_stat, 1);
    wt_stat = zeros(T_stat, 1);
    rt_stat = zeros(T_stat, 1);

    for t = 1:T_stat
        it_stat(t) = kt_stat(t + 1) - (1 - delta)*kt_stat(t);
        yt_stat(t) = exp(theta(shocks_stat(t)))*A*kt_stat(t)^alpha;
        ct_stat(t) = yt_stat(t) - it_stat(t);
        wt_stat(t) = (1 - alpha)*yt_stat(t);
        rt_stat(t) = alpha*yt_stat(t)/kt_stat(t);
    end

    kt_std(m) = std(log(kt_stat(100:T_stat)));
    yt_std(m) = std(log(yt_stat(100:T_stat)));
    ct_std(m) = std(log(ct_stat(100:T_stat)));
    it_std(m) = std(log(it_stat(100:T_stat)));
    wt_std(m) = std(log(wt_stat(100:T_stat)));
    rt_std(m) = std(log(rt_stat(100:T_stat)));

    kt_yt_correl(m) = corr(kt_stat(100:T_stat), yt_stat(100:T_stat));
    ct_yt_correl(m) = corr(ct_stat(100:T_stat), yt_stat(100:T_stat));
    it_yt_correl(m) = corr(it_stat(100:T_stat), yt_stat(100:T_stat));
    wt_yt_correl(m) = corr(wt_stat(100:T_stat), yt_stat(100:T_stat));
    rt_yt_correl(m) = corr(rt_stat(100:T_stat), yt_stat(100:T_stat));
end

disp(' ')
disp('Desviacion estandar del logaritmo: ')
disp(' ')

disp(['- Producto        = ', num2str(mean(yt_std))])
disp(['- Consumo         = ', num2str(mean(ct_std))])
disp(['- Inversión       = ', num2str(mean(it_std))])
disp(['- Capital         = ', num2str(mean(kt_std))])
disp(['- Salario real    = ', num2str(mean(wt_std))])
disp(['- Tasa de interés = ', num2str(mean(rt_std))])

disp(' ')
disp('Correlacion del producto y: ')
disp(' ')
disp(['- Consumo         = ', num2str(mean(ct_yt_correl))])
disp(['- Inversión       = ', num2str(mean(it_yt_correl))])
disp(['- Capital         = ', num2str(mean(kt_yt_correl))])
disp(['- Salario Real    = ', num2str(mean(wt_yt_correl))])
disp(['- Tasa de interés = ', num2str(mean(rt_yt_correl))])
disp(' ')

%%%%%%%%%%%%%%%    Inciso (c)    %%%%%%%%%%%%%%%

% Nota: la realización debe ser de epsilon, no de theta.
% Nota: viene del inciso 3 del laboratorio 2.

T         = 50;
theta     = [0 0.5 1];
shocks    = ones(T, 1)*2;
shocks(1) = 3;
kt        = zeros(T + 1, 1);
ind       = zeros(T, 1);

[aux, ind(1)] = min(k < k0);
kt(1) = k(ind(1));

for t = 1:T
   ind(t + 1) = G(ind(t), shocks(t));
   kt(t + 1)  = k(ind(t + 1));
end

it = zeros(T, 1);
ct = zeros(T, 1);
yt = zeros(T, 1);
wt = zeros(T, 1);
rt = zeros(T, 1);

for t = 1:T
    it(t) = kt(t + 1) - (1 - delta)*kt(t);
    yt(t) = exp(theta(shocks(t)))*A*kt(t)^alpha;
    ct(t) = yt(t) - it(t);
    wt(t) = exp(theta(shocks(t)))*A*(1 - alpha)*kt(t)^(alpha);
    rt(t) = exp(theta(shocks(t)))*A*(alpha)*kt(t)^(alpha - 1);
end

figure(3)

subplot(3, 2, 1);
plot(1:T, theta(shocks));
title('Choque Tecnológico');

subplot(3, 2, 2);
plot(1:T + 1, kt);
title('Capital');

subplot(3, 2, 3);
plot(1:T, yt);
title('Producto');

subplot(3, 2, 4);
plot(1:T, ct);
title('Consumo');

subplot(3, 2, 5);
plot(1:T, wt);
title('Salario Real');

subplot(3, 2, 6);
plot(1:T, rt);
title('Tasa de interes');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% De aquí para abajo ya no sirve.

%% Inciso 3

estado_choque = [1 2 3];
Pi            = [0.5 0.3 0.2;
                 0.1 0.7 0.2;
                 0.1 0.4 0.5];
Pi0 = [0 1 0];
T   = 50;
[X_path, states] = markov(estado_choque, Pi, Pi0, T);
shocks = X_path';

k0    = exp(0.5)*kss;
kt    = zeros(T + 1, 1);
ind   = zeros(T, 1);

[aux, ind(1)] = min(k < k0);
kt(1) = k(ind(1));

for t = 1:T
   ind(t + 1) = G(ind(t), shocks(t));
   kt(t + 1)  = k(ind(t + 1));
end

it = zeros(T, 1);
ct = zeros(T, 1);
yt = zeros(T, 1);
wt = zeros(T, 1);
rt = zeros(T, 1);

for t = 1:T
    it(t) = kt(t + 1) - (1 - delta)*kt(t);
    yt(t) = exp(theta(shocks(t)))*A*kt(t)^alpha;
    ct(t) = yt(t) - it(t);
    wt(t) = exp(theta(shocks(t)))*A*(1 - alpha)*kt(t)^(alpha);
    rt(t) = exp(theta(shocks(t)))*A*(alpha)*kt(t)^(alpha - 1);
end

figure(2)

subplot(3, 2, 1);
plot(1:T, theta(shocks));
title('Choque tecnologico');

subplot(3, 2, 2);
plot(1:(T + 1), kt);
title('Capital');

subplot(3, 2, 3);
plot(1:T, yt);
title('Producto');

subplot(3, 2, 4);
plot(1:T,ct);
title('Consumo');

subplot(3, 2, 5);
plot(1:T,wt);
title('Salario real');

subplot(3, 2, 6);
plot(1:T,rt);
title('Tasa de interes');
