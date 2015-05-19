function [ num ] = piso( g, p )

% num es el índice del trabajo óptimo
% La función piso tiene como argumentos la posición óptima para cada 
% combinación de estados y el tamaño de la malla del capital. 
% La posición g es un número natural entre 1 y p*r, este número
% tenemos que traducirlo a la posición óptima para el trabajo,
% es decir, transformarlo en un número natural entre 1 y r.
% Si el módulo de g y p es cero, entonces la función nos regresa el
% máximo entero menor que el cociente g/p. De lo contrario, le pedimos
% que nos recorra la posición del trabajo uno hacia adelante.

if mod(g,p)==0
    num=floor(g/p);
else
    num=floor(g/p)+1;
end

end

