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

A     = 1;      % Productividad
alpha = 0.40;   %
beta  = 0.987;  % Factor de descuento
delta = 0.012;  % 
gamma = 0.64;   % 

%% Parámetros de la solución

max_iter = 1000;  % Iteraciones máximas
tol      = 1e-2;  % Tolerancia
p        = 500;   % Malla

%%%%%%%%%%%%%%%    Inciso (a)    %%%%%%%%%%%%%%%

%
% Aproximación con el método de Tauchen (1986)
%

% Cinco estados para el shock
theta = [0.0231 -0.0115 0 0.0115 0.0231]';
q     = length(theta);

% Matriz de transición
Pi = [0.9727 0.0273      0      0      0;
      0.0041 0.9806 0.0153      0      0;
           0 0.0082 0.9837 0.0082      0;
           0      0 0.0153 0.9806 0.0041;
           0      0      0 0.0273 0.9727];

% Trabajo del estado estacionario
lss  = 1/(((gamma + (1 - alpha)*(1 - gamma))*(1 - beta*(1 - delta)) - ...
           alpha*beta*delta*gamma) / ...
          ((1 - alpha)*(1 - gamma)*(1 - beta*(1 - delta))));

% Capital del estado estacionario
kss = ((A*beta*alpha) / (1 - (1 - delta)*beta))^(1/(1 - alpha));

k     = linspace(2/3*kss, 1.5*kss, p);

% k0        = 2/3*kss;
k0        = zeros(1, p);
k1        = zeros(1, p);
% k(1)      = 0.000001;
l0        = zeros(q, p);
k1(p)     = k(p);
l_inicial = 0.8*lss;

for i = 1:p-1
    k1(i) = k(i+1);
end

for i = 1:p
    k0(i) = k(i);
end

%% TODO DE AQUI (Error: me da valores imaginarios)
options = optimset('Display','off');

for i = 1:p
    for j = 1:q
        z = theta(j);
        % Condición de eficiencia laboral
        l0(j,i) = fsolve(@(l0)cond_eficiencia_laboral(l0, ...
                                                 k0(i), ...
                                                 k1(i), ...
                                                 A, ...
                                                 alpha, ...
                                                 delta, ...
                                                 gamma, ...
                                                 z), l_inicial, options);
    end
end
valor = zeros(p*q, p, q);
%% HASTA AQUI

% Matriz valor: k(t) (capital)
for i = 1:p
    for j = 1:q
        valor((i - 1)*q + j, :, 1) = k(i)*ones(1, p);
    end
end

% Matriz valor: z(t) (shock)
vector = theta;
for i = 1:p-1
    vector = [vector; theta];
end
valor(:, :, 2) = vector*ones(1, p);

% Matriz valor: k(t+1) (capital)
valor(:, :, 3) = ones(q*p, 1)*k;

%% TODO posible fuente de error
% Problema: mucha volatilidad
% Queremos mismo approach que con theta

% Matriz valor: l(t) (trabajo)
locol = zeros(q*p, 1);
for i = 1:p
    locol(i)       = l0(1, i);
    locol(i + p)   = l0(2, i);
    locol(2*p + i) = l0(3, i);
    locol(3*p + i) = l0(4, i);
    locol(4*p + i) = l0(5, i);
end
val(:, :, 4) = locol*ones(1, p);
%% HASTA AQUI

M = zeros(p*q, p);

%% TODO posible fuente de error
% Generar matriz M y eliminar elementos imposibles
M = log(max(exp(val(:,:,2)).*A.* ...
            val(:,:,1).^alpha.* ...
            val(:,:,4).^(1 - alpha) + ...
            (1 - delta).*(val(:,:,1)) - val(:,:,3), ...
            1e-10));
%% HASTA AQUI

%
% Ejecución del algoritmo
%

E = eye(q);
for i = 1:p-1
    E = [E; eye(q)];
end

V0         = zeros(p*q, 1);
finalizado = 0;
iter       = 1;

while finalizado == 0 && iter < max_iter
    [V1, G1] = max((M + beta*(E*(Pi*reshape(V0, q, p))))');
    V1 = V1';
    G1 = G1';
    if norm(V0 - V1) < tol
        finalizado = 1;
    end
    disp(['Iter = ', num2str(iter), ...
          ' Error = ', num2str(norm(V0 - V1))])
    V0 = V1;
    iter = iter + 1;
end

%
% Extracción de resultados
%
G  = zeros(p, q);
V  = zeros(p, q);

for i = 1:p
    for j = 1:q
        G(i, j) = G1((i - 1)*q + j);
        V(i, j) = V1((i - 1)*q + j);
    end
end

% Valores iniciales
for j = 1:q
    G(1, j) = 1;
    V(1, j) = V(2, j);
end

c = zeros(p, q);

%% TODO Posible fuente de error
for i = 1:p
    for j = 1:q
        c(i, 1) = exp(theta(j))*A*k(i)^alpha*l0(1, j)^(1 - alpha) - ...
                  k(G(i, j)) + (1 - delta)*k(i);
    end
end
%% HASTA AQUI

%
% Gráficas
%

% figure(1);
% subplot(2, 1, 1);
% plot(k, k(G(:, 1)), 'b', ...
%      k, k(G(:, 2)), 'r', ...
%      k, k(G(:, 3)), 'g');
%      legend('theta = 0', 'theta = 0.5', 'theta = 1');
%      title('Capital optimo');

% subplot(2, 1, 2);
% plot(k, c_1, 'b', ...
%      k, c_2, 'r', ...
%      k, c_3, 'g');
%      legend('theta = 0', 'theta = 0.5', 'theta = 1');
%      title('Consumo optimo');
     
leyenda = legend('theta = -0.0231', ...
                 'theta = -0.0115', ...
                 'theta = 0', ...
                 'theta = 0.0115', ...
                 'theta = 0.0231');

% (a) Reglas de politica optima para el capital
figure(1)
hold all
plot(k(G(:, 1)))
plot(k(G(:, 2)))
plot(k(G(:, 3)))
plot(k(G(:, 4)))
plot(k(G(:, 5)))
set(leyenda, 'Location', 'SouthEast');
xlabel({'Periodos'});
ylabel({'k'});
axis([0 p 24 45]);
title('Reglas de politica optima para el capital');

% (c) Consumo
figure(2)
hold all
plot(c(:, 1))
plot(c(:, 2))
plot(c(:, 3))
plot(c(:, 4))
plot(c(:, 5))
set(leyenda,'Location','SouthEast');
xlabel({'Periodos'});
ylabel({'Consumo'});
axis([0 p 0.6 2.3]);
title('Trayectoria del consumo para distintos estados');
disp('Las graficas de las reglas de politica optima para k(t+1) y c(t) son:');
     
%%%%%%%%%%%%%%%    Inciso (b)    %%%%%%%%%%%%%%%

% Nota: Ajustar; viene del inciso 5 del laboratorio 2.


%%%%%%%%%%%%%%%    Inciso (c)    %%%%%%%%%%%%%%%

% Nota: Ajustar, viene del inciso 3 del laboratorio 2


%%%%%%%%%%%%%%%    Inciso (d)    %%%%%%%%%%%%%%%

%% TODO Revisar de aquí hasta el final

% Este archivo calcula las variables en estado estacionario;
% se usa para calibrar un nivel de gama y de impuestos que den como resultado
% las horas trabajadas en Mexico. Se supone un impuesto a los ingresos
% laborales, regresado a los agentes como transferencia lumpsum.
clear
clc

% Parametros del modelo
A     = 1;
alpha = 0.40;
beta  = 0.987;
delta = 0.012;
gama  = 0.64;
tauy  = 0;

% Calculo de variables en estado estacionario
lss = 1/(((gama+(1-alpha)*(1-tauy)*(1-gama))*(1-beta*(1-delta))-alpha*beta*delta*gama)/((1-alpha)*(1-gama)*(1-tauy)*(1-beta*(1-delta))));
kss = lss*((A*beta*alpha)/((1-beta*(1-delta))))^(1/(1-alpha));
yss = A*(kss^alpha)*(lss^(1-alpha));
iss = delta*kss;
css = yss-iss;
ratiokl = kss/lss;

% Ahora resolvemos la condicion de eficiencia del mercado laboral para
% obtener el valor de gama que resulta en lss = 0.45 (OECD, 2013), manteniendo k/l
% constante:
lss1 = 0.45;
kss1 = lss1*((A*beta*alpha)/((1-beta*(1-delta))))^(1/(1-alpha));
gamaprima = fsolve(@(gamaprima)eficienciamdol(gamaprima,tauy,lss1,kss1,A,alpha,delta),gama);

% Calculo de variables en estado estacionario con nuevo valor de gama
yss1 = A*(kss1^alpha)*(lss1^(1-alpha));
iss1 = delta*kss1;
css1 = yss1-iss1;
ratiokl1 = kss1/lss1;

% Ahora resolvemos la condicion de eficiencia del mercado laboral para
% obtener el valor de tauy que resulta en lss = 0.45, manteniendo k/l constante:
lss2 = lss1;
kss2 = kss1;
tauy1 = fsolve(@(tauy1)eficienciamdol1(tauy1,lss2,kss2,A,alpha,delta,gama),tauy);

% Calculo de variables en estado estacionario con nuevo valor de tauy
yss2 = A*(kss2^alpha)*(lss2^(1-alpha));
iss2 = delta*kss1;
css2 = yss2-iss2;
ratiokl2 = kss1/lss1;

% El hayazgo es consistente con Prescott 2004, a continuacion una curva que
% describe las ofertas de trabajo a distintas tasas impositivas;
grid = 100;
lssg = zeros(grid,1);
taus = linspace(-1.4,0.3,grid);

for i = 1:grid
    lssg(i) = 1/(((gama+(1-alpha)*(1-taus(i))*(1-gama))*(1-beta*(1-delta))-alpha*beta*delta*gama)/((1-alpha)*(1-gama)*(1-taus(i))*(1-beta*(1-delta))));
end

plot(taus,lssg);
axis([-1.4 0.3 0.226 0.5])
