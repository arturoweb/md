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

% alpha      = 0.40;   %
% beta       = 0.987;  % Factor de descuento
% delta      = 0.012;  % 

% TODO (otn) ¿dónde se deben usar estos parámetros?
% gamma      = 0.64;
% rho        = 0.95;
% sigma_eps  = 0.0024;
% sigma_zeta = 0.007;

%% Parámetros de la solución

max_iter = 1000;  % Iteraciones máximas
tol      = 1e-3;  % Tolerancia
A        = 1;     % Productividad
p        = 500;   % Malla

%%%%%%%%%%%%%%%    Inciso (a)    %%%%%%%%%%%%%%%

% Cuando esté listo para probarse el algoritmo
% debemos usar estos parámetros:

% Cinco estados para el shock
% theta = [0.0231 -0.0115 0.0115 0.0231]';
% q     = length(theta);

% Matriz de transición
% Pi = [0.9727 0.0273      0      0      0;
%       0.0041 0.9806 0.0153      0      0;
%            0 0.0082 0.9837 0.0082      0;
%            0      0 0.0153 0.9806 0.0041;
%            0      0      0 0.0273 0.9727];

% Con estos parámetros funciona bien:
alpha = 0.35;
beta  = 0.99;
delta = 0.06;
A     = 10;
theta = [0 0.5 1]';
q     = length(theta);
Pi    = [0.5 0.3 0.2;
         0.1 0.7 0.2;
         0.1 0.4 0.5];

% Capital de estado estacionario simple
lss  = 1/(((gamma + (1 - alpha)*(1 - gamma))*(1 - beta*(1 - delta)) - ...
           alpha*beta*delta*gamma) / ...
          ((1 - alpha)*(1 - gamma)*(1 - beta*(1 - delta))));
kss = ((A*beta*alpha) / (1 - (1 - delta)*beta))^(1/(1 - alpha));

%% TODO (otn) DIFERENTE DE AQUI

k0  = kss; % Estado estacionario vs malla

% TODO (otn) ¿tenemos que ajustar esto?
% k0  = (2/3)*kss;

% TODO (otn) ¿Por qué este tamaño y no otro?
% k = linspace(0, 1.2*exp(1)*kss, p);
k = linspace(0.8, 1.5*kss, p);

k(1) = 0.000001;

%% ESTO ES IMPORTANTE

% TODO (otn) Alejandro coloca aquí condiciones de eficiencia laboral

% for i = 1:p
%     for j = 1:q
%         z = theta(j);
%         l0(j,i) = fsolve(@(l0)seficiencialaboral(l0, k0(i), k1(i), A, alpha, ...
%                                                  delta, gamma, z), lini);
%     end
% end


%% HASTA AQUI

valor = zeros(p*q, p, q);

% Matriz valor, k(t)
for i = 1:p
    for j = 1:q
        valor((i - 1)*q + j, :, 1) = k(i)*ones(1, p);
    end
end

%% CONCILIAR DE AQUI

% Matriz valor, z(t) (vector theta)
% TODO (otn) No sé qué tiene Alejandro en 2
valor(:, :, 2) = ones(p*q, 1)*k;

% Matriz valor, k(t+1)
vector = theta;
for i = 1:p-1
    vector = [vector; theta];
end

% TODO (otn) alguien tiene las construcciones mal

% TODO (otn) Alejandro tiene k(t+1) en 3
valor(:, :, 3) = vector*ones(1, p);

% TODO (otn) Matriz valor, l(t) (en 4 ¿Por qué esta construcción?)
% locol = zeros(q*grid,1);
% for i = 1:grid
%     locol(i) = l0(1,i);
%     locol(i + grid)   = l0(2, i);
%     locol(2*grid + i) = l0(3, i);
%     locol(3*grid + i) = l0(4, i);
%     locol(4*grid + i) = l0(5, i);
% end
% val(:, :, 4) = locol*ones(1, grid);

% Matriz M
M = zeros(p*q, p);

% TODO (otn) tenemos esto diferente
M = log(max(exp(valor(:, :, 3)).*A.*valor(:, :,1).^(alpha) - ...
            valor(:, :, 2) + (1 - delta).*valor(:, :, 1), ...
            0));

%% HASTA AQUI

% Identidad repetida
E = eye(q);
for i = 1:p-1
    E = [E; eye(q)];
end

%
% Ejecución del algoritmo
%
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
    disp(['Iter = ', num2str(iter), ...
          ' Error= ', num2str(norm(V0 - V1))])
    V0 = V1;
    iter = iter + 1;
end

%
% Extracción de resultados
%
G1 = G;
G  = zeros(p, q);
V  = zeros(p, q);

for i = 1:p
    for j = 1:q
        G(i, 1) = G1((i - 1)*q + j);
        V(i, 3) = V1((i - 1)*q + j);
    end
end

% Valores iniciales
for j = 1:q
    G(1, j) = 1;
    V(1, j) = V(2, j);
end

% TODO (otn) modificar al método de Alejandro

for i = 1:p
    c_1(i) = exp(theta(1))*A*k(i)^alpha - k(G(i, 1)) + (1 - delta)*k(i);
    c_2(i) = exp(theta(2))*A*k(i)^alpha - k(G(i, 2)) + (1 - delta)*k(i);
    c_3(i) = exp(theta(3))*A*k(i)^alpha - k(G(i, 3)) + (1 - delta)*k(i);
end

%
% Gráficas
%

% TODO (otn) hacer las gráficas como Alejandro (más informativas)

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

% Nota: Ajustar; viene del inciso 5 del laboratorio 2.


%%%%%%%%%%%%%%%    Inciso (c)    %%%%%%%%%%%%%%%

% Nota: Ajustar, viene del inciso 3 del laboratorio 2


%%%%%%%%%%%%%%%    Inciso (d)    %%%%%%%%%%%%%%%

% Nota: Alejandro que pongamos un impuesto y me 
%       parece buena idea.