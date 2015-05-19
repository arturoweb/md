%   119399 David Alejandro Martell Juárez
%   118123 José Luis Cruz Álvarez
%   114370 Montserrat Ávila Acosta
%   Proyecto Final
%   Macroeconomía Dinámica 1
%   Profesor Carlos Urrutia
%  
%Modelo de ciclos económicos reales con decisiones de ocio-consumo.
%El Planificador Social resuelve la ecuación de Bellman: 
%       v(k,z) =  max {(1-gamma)*ln(exp(z)*k^alpha*l^(1-alpha)+(1-delta)*k-k'))+gamma*ln(1-l) + beta*E(v(k',z'))}
%                 l,k'
%       s.a 0<=k'<=exp(z)*k^alpha*l^(1-alpha)+(1-delta)*k
%           0<=l<=1
  
%
close 
clear all
clc

% Parámetros del modelo

alpha = 0.4;
gamma = 0.64;
beta  = 0.987;
delta = 0.012;
A=1;
rho=.95;
z = [-0.0231 -.0115 0 0.0115 0.0231];
expz = exp(z);
pi    = [0.9727 0.0273 0 0 0;0.0041 0.9806 0.0153 0 0;0 0.0082 0.9837 .0082 0; 0 0 .0153 .9806 .0041; 0 0 0 .0273 .9727];


% Otros parámetros requeridos en el programa

maxit = 1000; %Número máximo de iteraciones permitidas
p     = 170;  %Dimensión de la malla para valores del capital
crit  = 1e-2; %Criterio de tolerancia para la convergencia
q     = 5;    %Número de posibles estados del shock 
r=150;         %Dimensión de la malla para valores del trabajo

% Obtenemos el capital de estado estacionario (kss) para una economía
% sin incertidumbre (z=0) con el procedimiento especificado en el trabajo
% escrito.

phi = (1-beta)/beta;
gamma1 = gamma/(1-gamma);
sss = (alpha*delta)/(phi+delta);
lss = (1-alpha)/(1-alpha + gamma1*(1-sss));
kss   = ((alpha/(phi+delta))^(1/(1-alpha)))*lss;

% Se definen las mallas para el capital y el trabajo
k     = linspace(0.7*kss,1.2*kss,p);
l = linspace(0.25,0.35,r);

% Construimos la matriz M como se especifica en las notas de clase
% Se define una matriz de ceros en R4 para los posibles valores de (k,k',z,l)

value = zeros(p*q,p*r,4);
%Generamos una matriz con los posibles valores para el capital hoy (k)
F=k(1)*ones(q,p*r);
for i=2:p
    F = [F;k(i)*ones(q,p*r)];
end
value(:,:,1) = F;

%Generamos una matriz con los posibles valores para el capital mañana (k')
matriz=ones(p*q,1)*k;
for i=1:r-1
    matriz=[matriz ones(p*q,1)*k];
end
value(:,:,2)=matriz;

%Generamos una matriz con los posibles valores de los choques tecnológicos
%(z)
vece=z'*ones(1,p*r);
Q=z'*ones(1,p*r);
for i=1:p-1
    Q=[Q; vece];
end
value(:,:,3)=Q;

%Generamos una matriz con los posibles valores del trabajo (l)
vector2 = ones(p*q,p)*l(1);
for i=2:r
    vector2=[vector2 ones(p*q,p)*l(i)];
end
value(:,:,4)=vector2;

%Se llena la matriz M, con el valor de la función de utilidad para cada
%posible combinación (k,k',z,l)
M = (1-gamma)*log(max(exp(value(:,:,3)).*value(:,:,1).^(alpha).*value(:,:,4).^(1-alpha)+(1-delta).*value(:,:,1)-value(:,:,2),1e-320))+(gamma)*log(max(ones(p*q,p*r)-value(:,:,4),1e-320));
   
% Construimos la matriz auxiliar E  
    
I   = eye(q);
E   = I;
for i=1:p-1
    E=[E;I];
end


% Inicializamos la función valor (V0=0)
  V0 = zeros(p*q,1);

%Iteración de la función valor 
terminate = 0;
it = 1;
H = E*(pi*reshape(V0,q,p)); 

while (terminate==0 && it<maxit)
       
    H = E*(pi*reshape(V0,q,p));
    for i=1:r-1  
        H=[H E*(pi*reshape(V0,q,p))]; %Se construye la matriz H para 
                                      %ajustar las dimensiones 
    end  
	[V1,G] = max((M+beta*H)');
	V1     = V1';
	G      = G';
	if norm(V0-V1)<crit
		terminate=1;
	end
	disp(['Iteration= ',num2str(it),'  Error= ',num2str(norm(V0-V1))])
	V0     = V1;
	it     = it+1;
end

%Se obtienen las reglas de política óptima y la función valor

G1=reshape(G,q,p);
V = reshape(V1,q,p);

%Se obtienen las trayectorias óptimas para el capital, trabajo y consumo,
%dado un estado inicial (z,k)
%Ver script file para piso.m y resid.m
for i=1:p
    tra1=piso(G1(1,i),p);
    tra2=piso(G1(2,i),p);
    tra3=piso(G1(3,i),p);
    tra4=piso(G1(4,i),p);
    tra5=piso(G1(5,i),p);
    
    capi1=resid(G1(1,i),p);
    capi2=resid(G1(2,i),p);
    capi3=resid(G1(3,i),p);
    capi4=resid(G1(4,i),p);
    capi5=resid(G1(5,i),p);

    k_1(i) = k(capi1);
    k_2(i) = k(capi2);
    k_3(i) = k(capi3);
    k_4(i) = k(capi4);
    k_5(i) = k(capi5);

    l_1(i) = l(tra1);
    l_2(i) = l(tra2);
    l_3(i) = l(tra3);
    l_4(i) = l(tra4);
    l_5(i) = l(tra5);
    
    c_1(i)  = exp(z(1))*k(i)^alpha*l_1(i)^(1-alpha)+(1-delta)*k(i)-k_1(i);
    c_2(i)  = exp(z(2))*k(i)^alpha*l_2(i)^(1-alpha)+(1-delta)*k(i)-k_2(i);
    c_3(i)  = exp(z(3))*k(i)^alpha*l_3(i)^(1-alpha)+(1-delta)*k(i)-k_3(i);
    c_4(i)  = exp(z(4))*k(i)^alpha*l_4(i)^(1-alpha)+(1-delta)*k(i)-k_4(i);
    c_5(i)  = exp(z(5))*k(i)^alpha*l_5(i)^(1-alpha)+(1-delta)*k(i)-k_5(i);
end
   
%Se grafican las trayectorias óptimas de capital, función valor, consumo y
%trabajo para cada posible valor de los estados (z,k)

figure(1)
subplot(4,1,1)
plot(k,k_1,'b',k,k_2,'r',k,k_3,'g',k,k_4,'c',k,k_5,'m'), legend('z=-0.0231','z=-0.0115','z=0','z=0.0115','z=0.0231','Location','EastOutside'), title('Funciones de Política Óptima (Capital)')

subplot(4,1,2)
plot(k,V(1,:),'b',k,V(2,:),'r',k,V(3,:),'g',k,V(4,:),'c',k,V(5,:),'m'), legend('z=-0.0231','z=-0.0115','z=0','z=0.0115','z=0.0231','Location','EastOutside'), title('Función Valor')

subplot(4,1,3)
plot(k,c_1,'b',k,c_2,'r',k,c_3,'g',k,c_4,'c',k,c_5,'m'), legend('z=-0.0231','z=-0.0115','z=0','z=0.0115','z=0.0231','Location','EastOutside'), title('Consumo Óptimo')

subplot(4,1,4)
plot(k,l_1,'b',k,l_2,'r',k,l_3,'g',k,l_4,'c',k,l_5,'m'), legend('z=-0.0231','z=-0.0115','z=0','z=0.0115','z=0.0231','Location','EastOutside'), title('Trabajo Óptimo')

% Se simula una historia de shocks

T = 20; % Número de periodos deseados


% Se genera el vector de shocks AR(1)

shocks = ones(T,1);  % Se genera un vector para guardar los posibles 
                       % valores de los shocks.
shocks = shocks*3;     % Se da un valor arbitrario para el shock, el cual 
                       % puede cambiar voluntariamente.
shocks(10,1) = 4;                       

% Se simulan las trayectorias óptimas para el capital (k_t+1), trabajo
% (lt), inversión (it), consumo (ct) y producción (yt)

zt = zeros(T,1);
kt    = zeros(T+1,1);
posk   = zeros(T,1);
lt   = zeros(T,1);
posl   = zeros(T,1);
[aux,posk(1)] = min(k<kss); % Se parte del capital en la malla cuyo valor 
                            % es más cercano al del estado estacionario.
kt(1) = k(posk(1));

for t = 1:T
   posk(t+1) = resid(G1(shocks(t),posk(t)),p);
   kt(t+1)  = k(posk(t+1));
   posl(t) = piso(G1(shocks(t),posk(t)),p);
   lt(t) = l(posl(t));
   zt(t)=z(shocks(t));
end

it = zeros(T,1);
ct = zeros(T,1);
yt = zeros(T,1);

for t=1:T
    it(t) = kt(t+1)-(1-delta)*kt(t);
    yt(t) = expz(shocks(t))*kt(t)^alpha*lt(t)^(1-alpha);
    ct(t) = yt(t)-it(t);
end

% Se grafican las series de tiempo para las simulaciones anteriores

figure(2)
subplot(3,2,1), plot(1:T,zt), title('Choque Tecnológico')
subplot(3,2,2), plot(1:T+1,kt), title('Capital')
subplot(3,2,3), plot(1:T,yt), title('Producción')
subplot(3,2,4), plot(1:T,ct), title('Consumo')
subplot(3,2,5), plot(1:T,it), title('Inversión')
subplot(3,2,6), plot(1:T,lt), title('Trabajo')


% Calculamos los estadísticos

T_stat = 10000; %Número de periodos para la simulación

pi_acum = cumsum(pi',1)'; % Función de distribución.
                          % Otra manera de lograr esto es con 
                          % el siguiente comando: pi_acum1=cumsum(pi,2) 
                       
for m=1:100   % Número de iteraciones para calcular los estadísticos

shocks_stat = zeros(10000,1);
shocks_stat(1)=3;


for i=1:(T_stat)-1
    va = random('unif',0,1);
    for j=1:q
        for n=1:q
            if shocks_stat(i)==j
                if n==1
                    if va<=pi_acum(j,n)
                        shocks_stat(i+1)=n;
                    end
                end
                if n>1
                    if va<=pi_acum(j,n) && va>pi_acum(j,n-1)
                        shocks_stat(i+1)=n;
                    end
                end
            end
        end
    end
end

kt_stat    = zeros(T_stat+1,1);
posk_stat   = zeros(T_stat,1);
lt_stat   = zeros(T_stat,1);
posl_stat   = zeros(T_stat,1);
[aux, posk_stat(1)] = min(k<kss);
kt_stat(1) = k(posk_stat(1));

for t = 1:T_stat
   posk_stat(t+1) = resid(G1(shocks_stat(t),posk_stat(t)),p);
   kt_stat(t+1)  = k(posk_stat(t+1));
   posl_stat(t) = piso(G1(shocks_stat(t),posk_stat(t)),p);
   lt_stat(t) = l(posl_stat(t));
end

it_stat = zeros(T_stat,1);
ct_stat = zeros(T_stat,1);
yt_stat = zeros(T_stat,1);

for t=1:T_stat
    it_stat(t) = kt_stat(t+1)-(1-delta)*kt_stat(t);
    yt_stat(t) = expz(shocks_stat(t))*kt_stat(t)^alpha*lt_stat(t)^(1-alpha);
    ct_stat(t) = yt_stat(t)-it_stat(t);
end

% Se calculan las desviaciones estándar deseadas

kt_std(m) = std(log(kt_stat(1000:10000)));
yt_std(m) = std(log(yt_stat(1000:10000)));
ct_std(m) = std(log(ct_stat(1000:10000)));
it_std(m) = std(log(it_stat(1000:10000)));
lt_std(m) = std(log(lt_stat(1000:10000)));

% Se calculan las correlaciones deseadas
kt_yt_correl(m) = corr(kt_stat(1:T_stat),yt_stat);
ct_yt_correl(m)= corr(ct_stat,yt_stat);
it_yt_correl(m) = corr(it_stat,yt_stat);
lt_yt_correl(m) = corr(lt_stat,yt_stat);
m
end

disp(' ')
disp('Desviación Estándar del logaritmo de la serie: ')
disp(' ')

disp(['- Producción        = ',num2str(mean(yt_std))])
disp(['- Consumo   = ',num2str(mean(ct_std))])
disp(['- Inversión    = ',num2str(mean(it_std))])
disp(['- Capital       = ',num2str(mean(kt_std))])
disp(['- Trabajo       = ',num2str(mean(lt_std))])

disp(' ')
disp('Correlación de la producción con: ')
disp(' ')
disp(['- Consumo   = ',num2str(mean(ct_yt_correl))])
disp(['- Inversión    = ',num2str(mean(it_yt_correl))])
disp(['- Capital       = ',num2str(mean(kt_yt_correl))])
disp(['- Trabajo       = ',num2str(mean(lt_yt_correl))])
disp(' ')

figure(3)
subplot(3,2,1), plot(1:T_stat,z(shocks_stat)), title('Choque Tecnológico')
subplot(3,2,2), plot(1:T_stat+1,kt_stat), title('Capital')
subplot(3,2,3), plot(1:T_stat,yt_stat), title('Producción')
subplot(3,2,4), plot(1:T_stat,ct_stat), title('Consumo')
subplot(3,2,5), plot(1:T_stat,it_stat), title('Inversión')
subplot(3,2,6), plot(1:T_stat,lt_stat), title('Trabajo')