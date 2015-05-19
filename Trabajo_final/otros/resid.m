function [num] = resid(g,p)

% num es el �ndice de la posici�n del capital �ptimo
% La funci�n resid tiene como argumentos la posici�n �ptima para cada 
% combinaci�n de estados y el tama�o de la malla del capital.
% La posici�n g es un n�mero natural entre 1 y p*r, este n�mero
% tenemos que traducirlo a la posici�n �ptima del capital de ma�ana,
% es decir, transformarlo en un n�mero natural entre 1 y p. 
% Lo que este algoritmo realiza es obtener el modulo de la posici�n
% g entre el tama�o de la malla para el capital, si la divisi�n es
% exacta, entonces el capital �ptimo corresponde al �ltimo elemento de 
% la malla para un cierto valor del trabajo. Ahora bien, si el m�dulo es
% distinto de cero, entonces la posici�n tomar� este valor.

if (mod(g,p)~=0)
    num = mod(g,p);
else
    num = p;
end