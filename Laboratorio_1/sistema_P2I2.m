%
% Laboratorio 1
%
% Sistema de la pregunta 2, inciso 2 (P2I2).
%
% Omar Trejo Navarro, 119711
% Macroeconomia Dinamica I,
% Prof. Carlos Urrutia,
% ITAM, 2015
%
function [s] = sistema_P2I2(x)
    s = [ 2*x(1) - x(2) - exp(-x(1)),
          -x(1) - 2*x(2) - exp(-x(2)) ];
end
