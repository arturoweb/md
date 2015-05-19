function [num] = resid(g,p)

% num es el índice de la posición del capital óptimo
% La función resid tiene como argumentos la posición óptima para cada 
% combinación de estados y el tamaño de la malla del capital.
% La posición g es un número natural entre 1 y p*r, este número
% tenemos que traducirlo a la posición óptima del capital de mañana,
% es decir, transformarlo en un número natural entre 1 y p. 
% Lo que este algoritmo realiza es obtener el modulo de la posición
% g entre el tamaño de la malla para el capital, si la división es
% exacta, entonces el capital óptimo corresponde al último elemento de 
% la malla para un cierto valor del trabajo. Ahora bien, si el módulo es
% distinto de cero, entonces la posición tomará este valor.

if (mod(g,p)~=0)
    num = mod(g,p);
else
    num = p;
end