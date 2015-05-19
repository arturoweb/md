function [ num ] = piso( g, p )

% num es el �ndice del trabajo �ptimo
% La funci�n piso tiene como argumentos la posici�n �ptima para cada 
% combinaci�n de estados y el tama�o de la malla del capital. 
% La posici�n g es un n�mero natural entre 1 y p*r, este n�mero
% tenemos que traducirlo a la posici�n �ptima para el trabajo,
% es decir, transformarlo en un n�mero natural entre 1 y r.
% Si el m�dulo de g y p es cero, entonces la funci�n nos regresa el
% m�ximo entero menor que el cociente g/p. De lo contrario, le pedimos
% que nos recorra la posici�n del trabajo uno hacia adelante.

if mod(g,p)==0
    num=floor(g/p);
else
    num=floor(g/p)+1;
end

end

