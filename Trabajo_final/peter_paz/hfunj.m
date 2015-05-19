function F = hfunj(l,k,kp,z)

%=======================================================================
%Funcion que resuelve la condicion de primer orden del trabajo y consumo
%autor: PETER PAZ
%=======================================================================
% Para que la funcion sea usada necesita valores del capital,
% capital de mañana y choques.
%
A=1;
alpha = .4;
gamma =.64;
delta = 0.012;

y = exp(z) *A* k^alpha * l^(1-alpha);
I = kp - (1 - delta) * k;
c = y - I;

F = (1 - alpha) * (y / c)- (gamma/(1 - gamma)) * ((l) / (1 - l)) ;
%=======================================================================
%Fin de la funcion
%=======================================================================