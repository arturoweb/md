%% Resolviendo el ejercicio ii

% Parametros del modelo
alpha = .35;
beta  = .99;
delta = .06;
A     = 10;
theta = [0; 0.5; 1];

% Generando maya
maxit = 1000;
p     = 500;
crit  = 1e-4;
T     = 50;

kss = ((A*beta*alpha)/(1 - (1 - delta)*beta))^(1/1 - alpha));

% En el mejor de los casos, si el shock no nos afecta
% llegamos a este capital de estado estacionario.

k = linspace(0, 1.5*kss, p);

% Maya de k, tres matrices
val = zeros(p*q, p, 3);
for i = 1:p
    val(((i - 1)*q) + 1, :, 1) = k(i)*ones(1, p);
    val(((i - 1)*q) + 2, :, 1) = k(i)*ones(1, p);
    val(((i - 1)*q) + 3, :, 1) = k(i)*ones(1, p);
end

val(:, :, 2) = ones(p*q, 1)*k;
vector = theta;

for i=1:p-1
    vector=[vector;theta];
end

val(:, :, 3) = vector*ones(1, p);
M = zeros(p*q,p);
M = log(max(
        exp(val(:, :, 3)) .* A .* ...
        ((val(:, :, 1)).^alpha + ...
        (1 - delta).*val(:, :, 1) - val(:, :, 2), ...
        1e-40));

% Matriz auxiliar
I = eye(q);

% Generando E
for i = 1:p - 1
    E = [E,I];
end

% Iteración funcion de valor
v0 = zeros(p*q, 1);
fin = 0;
it = 1;

while (fin == p0 ite < maxit)
    [V1, G1] = max(M + beta*(E*(pi*reshape(V0, q, p))))');
    V1 = V1';
    G1 = G1';
    if norm(V0 - V1) < crit
        fin01;
    end
    V0 = V1;
    ite = ite+1;
end
