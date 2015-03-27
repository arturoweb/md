%
% Laboratorio 1
%
% Respuesta de la pregunta 4
%
% Omar Trejo Navarro, 119711
% Macroeconomia Dinamica I,
% Prof. Carlos Urrutia,
% ITAM, 2015
%

clear all
close all
clc

%% Parametros

alpha = 0.35;
beta  = 0.99;
delta = 0.06;
A     = 10;
tau   = 0;
maxit = 1000;
T     = 100;
kss   = ((A*beta*alpha)/(1 - (1 - delta)*beta))^(1/(1 - alpha));
k0    = (2/3)*kss;
malla = 500;
crit  = 1e-3;

%% Iteracion de la funcion valor

tic
[kt_vf, yt_vf, it_vf, ct_vf] = v(alpha, ...
                                 beta, ...
                                 delta, ...
                                 A, ...
                                 maxit, ...
                                 crit, ...
                                 T, ...
                                 kss, ...
                                 k0, ...
                                 malla);
toc

rt_vf = alpha.*A*kt_vf.^(alpha - 1); % Tasa de interés
wt_vf = (1-alpha).*A*kt_vf.^alpha; % Salarios

%% Gráficas

grafica_P3P4(kt_vf, yt_vf, it_vf, ct_vf, rt_vf, wt_vf, T);