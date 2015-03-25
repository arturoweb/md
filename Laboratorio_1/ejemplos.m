% Laboratorio 1 de Macro Din�mica I
% Laboratorista: Cristi�n Aguilera Arellano
% E-mail: ecointernacional2015@gmail.com

% Introducci�n

clc             % permite eliminar lo que est� en el command window
clear all       % permite borrar todo lo que se almacena en mi "Workspace"


    % Introducci�n

    % Generaci�n de Vectores
        %Generar un vector de 1x3
        A=[10 9 8];      % si ponemos ; al final la ejecuta sin mostrar en el command window.

        %Generar otro vector de 1x3
        B=[5 5 5];


        % Sumar ambos vectores
        D=A+B;

        % El m�ximo del vector D
        max(D);
        % El m�nimo del vector D
        min (D);
        % Sumar todos los elementos del vector B
        display('La suma de todos los elementos del vector B es:')
        sum(B);

        % Generando una matriz de 2x2
        Y=[1 2; 3 4];
        % Transponer la matriz de 2x2
        YT=Y';

        % Generaci�n de un vector con n�meros consecutivos del uno al diez
        Z=1:10;
        % Generaci�n de un vector con n�meros consecutivos del uno al diez
        % espaciando de dos en dos
        T=1:2:10;
        T=1:3:10;
        T=1994:2013; % El valor de T se va reescribiendo, lo que importa es el valor final.
        T(2);        % Muestra el segundo elemento del vector, en general T(n) muestra el n-�simo t�rmino del vector

        % Generando un vector de ceros
        K=zeros(10,1);
        K(2)=1;

    % Loops

        % For
        % Law motion of capital I(t)=K(t+1)-(1-delta)*K(t)

        % Primero inventemos una serie para la inversi�n
        % Supongamos que la inversi�n en el a�o 2000 es 1,  y que la tasa
        % de crecimiento de la inversi�n es 10%

        g=.1;    % Tasa de crecimiento de la inversi�n
        I(1)=1;

        for t=1:14

            I(t+1)=(1+g)*I(t);

        end

       K=zeros(15,1);
       K(1)=1;
       delta=.05;

       for t=1:14

           K(t+1)=I(t)+(1-delta)*K(t);

       end

    % Funciones
             % File -> New -> Function
             % Creamos una funci�n para la media y la varianza
             % Generamos un vector X con 10 elementos
             G=1:10;
             [ave,std]=stat(G);

    % Fsolve - resuelve ecuaciones no lineales.
            % Estructura: [x,fval]=fsolve(@myfun, x0)

