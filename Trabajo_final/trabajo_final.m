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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Del correo de Alejandro.. problemas:

% 1. Según yo, cuando k0<kss, la trayectoria hacia el ss está 
% caracterizada por lt decreciente. Me sale creciente. :´(

% Respueta: Según yo esto está bien. Si el capital empieza
% por debajo del de estado estacionario, entonces deberá ir
% aumentando, lo que a su vez hará más rentable el trabajo
% en el margen y provocará que suba el empleo.

% 2. Tanto en la descripción del lab 1 como en el m file de paz, le 
% dedican unas líneas a quitar las soluciones de lt que son mayores 
% a 1. A mí no me salen soluciones mayores a 1. Nunca. Pero estoy 
% seguro que hay algún error porque siempre la solución de lt está 
% por encima de lss.

% Respuesta: Comno l_t es proporcional, entonces si obtiene valor
% arriba de 1 está mal (y se debe eliminar, aunque no debería de
% ocurrir en primer lugar). Que no te salga l_t > 1 está bien,
% pero habrá que ver por qué no llega al de estado estacionario.

% 3. Chequé mi álgebra, para ver que no fuera un pedo de aritmética 
% con los parámetros, pero lss y kss ya cuadran con lo de paz. Tons 
% no es eso.

% Respuesta: No entendí esto.

% Archivos:
% t_final_6.m
% seficiencialaboral.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;
clc;

%% Parámetros
beta       = 0.987;
gamma      = 0.64;
alpha      = 0.40;
delta      = 0.012;
rho        = 0.95;
sigma_eps  = 0.0024;
sigma_zeta = 0.007;

max_iter   = 500;
grid       = 500;
tol        = 1e-3;
T          = 100;
A          = 1;     % ? A = 10
kss        = ((A*beta*alpha) / ...
              (1 - (1 - delta)*beta))^(1/(1 - alpha));
k0         = (2/3)*kss;
malla      = 500;

%% Inciso a

% Usando los siguientes valores de los parámetros (frecuencia 
% trimestral) resolver el problema del planificador social para 
% esta economía usando el método de iteración de la función de 
% valor y encontrar las reglas de decisión óptimas.

% NOTA: Pueden utilizar la siguiente aproximación discreta al 
% proceso AR(1), obtenida usando el método de Tauchen "Finite 
% State Markov Chain Approximations to Univariate and Vector 
% Autoregressionsion", Economic Letters, 20, 1986.
% Esta aproximación está caracterizada por cinco posibles estados 
% para el shock z_t:

% Cinco estados para el shock
theta = [0.0231 -0.0115 0.0115 0.0231];

% Matriz de transición
Pi = [0.9727 0.0273      0      0      0;
      0.0041 0.9806 0.0153      0      0;
           0 0.0082 0.9837 0.0082      0;
           0      0 0.0153 0.9806 0.0041;
           0      0      0 0.0273 0.9727];

tic
[kt_vf, yt_vf, it_vf, ct_vf] = iter_func_valor(alpha, ...
                                               beta, ...
                                               delta, ...
                                               A, ...
                                               max_iter, ...
                                               tol, ...
                                               T, ...
                                               kss, ...
                                               k0, ...
                                               malla);
toc

rt_vf = alpha.*A*kt_vf.^(alpha - 1);  % Tasa de interés
wt_vf = (1-alpha).*A*kt_vf.^alpha;    % Salarios

grafica_P3P4(kt_vf, yt_vf, it_vf, ct_vf, rt_vf, wt_vf, T);

%% Inciso b

% Usando una simulación larga (10; 000 perÌodos) parta aproximar 
% la distribución invariante, calcular la desviación estándar del 
% logaritmo de cada serie (en %) y su correlaciÛn con el producto. 
% Compare sus resultados con Cooley-Prescott (1995), Tabla 1.2.

% NOTA: Los resultados pueden no coincidir exactamente, puesto que 
% el ejercicio omite algunas caracterÌsticas del modelo de RBC
% (crecimiento exÛgeno de la tecnología y la población). Puesto 
% que la tendencia en el modelo es una constante, como primera 
% aproximación no es necesario Öltrar las series resultantes.

%% Inciso c

% Simular y graÖcar las funciones de impulso-respuesta para cada 
% series, es decir, la evoluciÛn de cada serie luego de una
% realización positiva de " (de una desviación est·ndar), 
% partiendo de un equilibrio estacionario.

%% Inciso d

% Cambiando alguno de los par·metros del modelo (que ustedes 
% consideren interesante), volver a ejecutar las partes (a) y 
% (b) y comparar los resultados de las tablas. Interpretar las 
% diferencias en resultados.