% Macro dinamica I
% Trabajo final

close all
clear all
clc

% Matriz de transicion
pi = [0.9727 0.0273 0 0 0; 0.041 0.9806 0.0153 0 0;0 0.0082 0.9837 0.0082 0;0 0 0.0153 0.9806 0.0041;0 0 0 0.0273 0.9727];

% Vector de shocks z(t)
theta = [-0.0231 -0.0115 0 0.0115 0.0231]';

% Parametros del modelo
A     = 1;    % maximo valor de z(t) para calcular edo. est.
alpha = 0.40;
beta  = 0.987;
delta = 0.012;
gamma = 0.64;

% Parametros de la solucion
miter = 500;
grid = 500;
tol = 1e-04;
q = length(theta);


% Malla de valores para el capital (ahora depende de lss):
lss  = 1/(((gamma+(1-alpha)*(1-gamma))*(1-beta*(1-delta))-alpha*beta*delta*gamma)/((1-alpha)*(1-gamma)*(1-beta*(1-delta))));
kss  = lss*((A*beta*alpha)/((1-beta*(1-delta))))^(1/(1-alpha));
kmin = 0.8*kss;
kmax = 1.5*kss;
k    = linspace(kmin,kmax,grid);

%Solucion de la condicion de eficiencia para obtener una malla de l(t)

%malla auxiliar de valores de k(t) y k(t+1) para fsolve:
k1 = zeros(1,grid);

for i = 1:grid-1
    k1(i) = k(i+1);
end

k1(grid) = k(grid);
k0 = zeros(1,grid);

for i = 1:grid
    k0(i) = k(i);
end

%solucion secuencial de la condicion de eficiencia en el mdo. laboral.
l0   = zeros(q,grid);
lini = 0.8*lss;

for i = 1:grid
    for j = 1:q
        z = theta(j);
        l0(j,i) = fsolve(@(l0)seficiencialaboral(l0, k0(i), k1(i), A, alpha, ...
                                                 delta, gamma, z), lini);
    end
end

% Generar la matriz val donde cada componente tiene cuatro elementos:
val = zeros(q*grid,grid,4);

% Primer elemento de cada componente en la matriz val: k(t)
for i = 1:grid
    val(((i-1)*q)+1,:,1) = k(i)*ones(1,grid);
    val(((i-1)*q)+2,:,1) = k(i)*ones(1,grid);
    val(((i-1)*q)+3,:,1) = k(i)*ones(1,grid);
    val(((i-1)*q)+4,:,1) = k(i)*ones(1,grid);
    val(((i-1)*q)+5,:,1) = k(i)*ones(1,grid);
end

% Segundo elemento de la matriz val: z(t) (elementos del vector theta)
vector = theta;

for i = 1:grid-1
    vector = [vector;theta];
end

val(:,:,2) = vector*ones(1,grid);

%Tercer elemento de cada componente en la matriz val: k(t+1)
val(:,:,3) = ones(q*grid,1)*k;

%Cuarto elemento de cada componente en la matriz val: l(t)
%pasar l0 a vector columna
locol = zeros(q*grid,1);
for i = 1:grid
    locol(i) = l0(1,i);
    locol(i+grid) = l0(2,i);
    locol(2*grid+i) = l0(3,i);
    locol(3*grid+i) = l0(4,i);
    locol(4*grid+i) = l0(5,i);
end;

val(:,:,4) = locol*ones(1,grid);

% Generar la matriz M y eliminar los elementos imposibles al usar max
M = zeros(q*grid,grid);
M = log(max((exp(val(:,:,2)).*(val(:,:,1)).^alpha.*(val(:,:,4)).^(1-alpha)+(1-delta).*(val(:,:,1))-val(:,:,3)),1e-40));

% Generar matriz auxiliar E
E = eye(q);
for i = 1:grid-1
    E = [E;eye(q)];
end

% Iterar la funcion de valor
V0 = zeros(q*grid,1);

fin = 0;
ite = 1;

while (fin =  = 0 && ite<miter)
    [V1,G1] = max((M+beta*(E*(pi*reshape(V0,q,grid))))');
    V1 = V1';
    G1 = G1';
    if norm(V1-V0)<tol
        fin = 1;
    end
    V0 = V1;
    ite = ite+1;
end

%Para obtener las reglas de politica optimas y la funcion de valor se construyen las siguientes matrices:
V = zeros(grid,q);
G = zeros(grid,q);

for i = 1:grid
    G(i,1) = G1((i-1)*q+1);
    G(i,2) = G1((i-1)*q+2);
    G(i,3) = G1((i-1)*q+3);
    G(i,4) = G1((i-1)*q+4);
    G(i,5) = G1((i-1)*q+5);

    V(i,1) = V1((i-1)*q+1);
    V(i,2) = V1((i-1)*q+2);
    V(i,3) = V1((i-1)*q+3);
    V(i,4) = V1((i-1)*q+4);
    V(i,5) = V1((i-1)*q+5);
end

% Ajuste del primer valor de G para los 3 estados en 1
G(1,1) = 1;
G(1,2) = 1;
G(1,3) = 1;
G(1,4) = 1;
G(1,5) = 1;

% Ajustamos del primer valor de V para los 3 estados igual al segundo
V(1,1) = V(2,1);
V(1,2) = V(2,2);
V(1,3) = V(2,3);
V(1,4) = V(2,2);
V(1,5) = V(2,3);

%Para obtener los valores optimos del consumo se construye la siguiente matriz:
c = zeros(grid,q);

for i = 1:grid
    c(i,1) = exp(theta(1))*k(i)^alpha*l0(1,i)^(1-alpha)-k(G(i,1))+(1-delta)*k(i);
    c(i,2) = exp(theta(2))*k(i)^alpha*l0(2,i)^(1-alpha)-k(G(i,2))+(1-delta)*k(i);
    c(i,3) = exp(theta(3))*k(i)^alpha*l0(3,i)^(1-alpha)-k(G(i,3))+(1-delta)*k(i);
    c(i,4) = exp(theta(4))*k(i)^alpha*l0(4,i)^(1-alpha)-k(G(i,4))+(1-delta)*k(i);
    c(i,5) = exp(theta(5))*k(i)^alpha*l0(5,i)^(1-alpha)-k(G(i,5))+(1-delta)*k(i);
end

%Graficos relevantes:

%(a) Reglas de politica optima para el capital
figure
hold all
%Capital en el estado 1
plot(k(G(:,1)))
%Capital en el estado 2
plot(k(G(:,2)))
%Capital en el estado 3
plot(k(G(:,3)))
%Capital en el estado 4
plot(k(G(:,4)))
%Capital en el estado 5
plot(k(G(:,5)))
leyenda = legend('theta = -0.0231','theta = -0.0115','theta = 0','theta = 0.0115','theta = 0.0231');
set(leyenda,'Location','SouthEast');
xlabel({'Periodos'});
ylabel({'k'});
axis([0 grid 24 45]);
title('Reglas de politica optima para el capital')

%(c) Consumo
figure
hold all
%Consumo en el estado 1
plot(c(:,1))
%Consumo en el estado 2
plot(c(:,2))
%Consumo en el estado 3
plot(c(:,3))
%Consumo en el estado 4
plot(c(:,4))
%Consumo en el estado 5
plot(c(:,5))
leyenda = legend('theta = -0.0231','theta = -0.0115','theta = 0','theta = 0.0115','theta = 0.0231');
set(leyenda,'Location','SouthEast');
xlabel({'Periodos'});
ylabel({'Consumo'});
axis([0 grid 0.6 2.3]);
title('Trayectoria del consumo para distintos estados')

disp('Las graficas de las reglas de politica optima para k(t+1) y c(t) son:')
