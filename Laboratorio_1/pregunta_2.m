%
% Laboratorio 1
%
% Respuesta de la pregunta 2
%
% Omar Trejo Navarro, 119711
% Macroeconomia Dinamica I,
% Prof. Carlos Urrutia,
% ITAM, 2015
%

clear all
close all
clc
options = optimset('Display', 'off');

%% Inciso 1

disp('Pregunta 2, inciso 1.');
x0 = [1 1]';
[x, fval] = fsolve(@(x) sistema_P2I1(x), x0, options);

str2 = sprintf(['El resultado es (x_1, x_2) = (%2.3f, %2.3f) '...
               'partiendo del punto inicial (x_1, x_2) = (%2.3f, %2.3f).'], ...
               x(1), x(2), x0(1), x0(2));

disp(str2);

%% Inciso 2

disp('Pregunta 2, inciso 2.');
x0 = [-5 -5]';
[x, fval] = fsolve(@(x) sistema_P2I2(x), x0, options);

str1 = sprintf(['El resultado es (x_1, x_2) = (%2.3f, %2.3f) ' ... 
                'partiendo del punto inicial (x_1, x_2) = (%2.3f, %2.3f).'], ...
                x(1), x(2), x0(1), x0(2));

disp(str1);

