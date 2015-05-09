%
% Laboratorio 2
%
% Omar Trejo Navarro, 119711
% Macroeconomia Dinamica I,
% Prof. Carlos Urrutia,
% ITAM, 2015
%

clear all
close all
clc

%% Inciso 1

P = [0.5 0.3 0.2;
     0.1 0.7 0.2;
     0.1 0.4 0.5];

A = P^1000;

%% Inciso 2

alpha = 0.35;
beta  = 0.99;
delta = 0.06;
A     = 10;
maxit = 1000;
p     = 500;
crit  = 1e-2; % Tolerancia
T     = 50;
q     = 3;    % 3 estados

theta = [0 0.5 1];
Pi    = [0.5 0.3 0.2;
         0.1 0.7 0.2;
         0.1 0.4 0.5];

% Estado estacionario sin incertidumbre
kss  = ((A*beta*alpha)/(1 - (1 - delta)*beta))^(1/(1 - alpha));
k0   = kss;
k    = linspace(0, 1.2*exp(1)*kss, p);
k(1) = 0.000001;

value = zeros(p*q,p,3);
for i=1:p
    value((i-1)*q+1,:,1) = k(i)*ones(1,p);
    value((i-1)*q+2,:,1) = k(i)*ones(1,p);
    value((i-1)*q+3,:,1) = k(i)*ones(1,p);
end

% Matriz de capitales por filas
value(:, :, 2) = ones(p*q, 1)*k;
vector = theta';
for i = 1:p-1
    vector = [vector; theta'];
end
value(:, :, 3) = vector*ones(1, p);
M      = zeros(p*q,p);
Mzeros = zeros(p*q,p);

%     F(Sij,Kl)=F(Ki,Zj,Kl)=u(exp(Zj)*f(Ki)+(1-delta)Ki-Kl)
M = log(max(exp(value(:,:,3)).*A.*value(:,:,1).^(alpha)-value(:,:,2)+...
    (1-delta).*value(:,:,1),0));

I = eye(q);
E = I;
for i = 1:p-1
    E = [E; I];
end

% Algoritmo
V0 = zeros(p*q, 1);
terminate = 0; % == 1 -> finalizado
it = 1;
while terminate == 0 && it < maxit
    [V1,G] = max((M + beta*(E*(Pi*reshape(V0, q, p))))');
    V1 = V1';
    G  = G';
    if norm(V0 - V1) < crit
        terminate = 1;
    end
    disp(['It = ', num2str(it), ' Error= ', num2str(norm(V0 - V1))])
    V0 = V1;
    it = it + 1;
end

G1 = G;
G = zeros(p,q);
V = zeros(p,q);

for i=1:p
    G(i,1) = G1((i-1)*q+1);
    G(i,2) = G1((i-1)*q+2);
    G(i,3) = G1((i-1)*q+3);
    V(i,1) = V1((i-1)*q+1);
    V(i,2) = V1((i-1)*q+2);
    V(i,3) = V1((i-1)*q+3);

end

G(1,1) = 1;
G(1,2) = 1;
G(1,3) = 1;
V(1,1) = V(2,1);
V(1,2) = V(2,2);
V(1,3) = V(2,3);

for i = 1:p
    c_1(i) = exp(theta(1))*A*k(i)^alpha-k(G(i, 1)) + (1 - delta)*k(i);
    c_2(i) = exp(theta(2))*A*k(i)^alpha-k(G(i, 2)) + (1 - delta)*k(i);
    c_3(i) = exp(theta(3))*A*k(i)^alpha-k(G(i, 3)) + (1 - delta)*k(i);
end

figure(1);

subplot(2, 1, 1);
plot(k, k(G(:, 1)), 'b', ...
     k, k(G(:, 2)), 'r', ...
     k, k(G(:, 3)), 'g'), ...
     legend('theta = 0', 'theta = 0.5', 'theta = 1'), ...
     title('Capital optimo');

subplot(2, 1, 2);
plot(k, c_1, 'b', ...
     k, c_2, 'r', ...
     k, c_3, 'g'), ...
     legend('theta = 0', 'theta = 0.5', 'theta = 1'), ...
     title('Consumo optimo');

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

%% Inciso 4

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

%% Inciso 5

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
