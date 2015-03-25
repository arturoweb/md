%
% Ejercicio 1
%

%
% Inciso 1
%

% Ecuacion 1
x0 = 1;
[x, f1] = fsolve(@funcion_P1I1E1, x0);
disp(x);

% Ecuacion 2
x0 = 1;
[x, f2] = fsolve(@funcion_P1I1E2, x0);
disp(x);

% Ecuacion 3
x0 = 1;
[x, f3] = fsolve(@funcion_P1I1E3, x0);
disp(x);

%
% Inciso 2
%
options = optimset('Display','off');

disp('Case 1');
x0 = 1.625;
[x, f4, flag] = fsolve(@funcion_P1I2, x0, options);

if flag == 1
    str = sprintf('El resultado es [%2.3f] partiendo del guess [x0=1.625]', x);
    disp(str);
else
    str = sprintf('El resultado no converge partiendo del guess [x0=1.625]');
    disp(str);
end

disp('Case 2');
x0 = 1.875;
[x2, fval5, flag] = fsolve(@funcion_P1I2, x0, options);

if flag == 1
    str = sprintf('El resultado es [%2.3f] partiendo del guess [x0=1.875]', x2);
    disp(str);
else
    str = sprintf('El resultado no converge partiendo del guess [x0=1.875]');
    disp(str);
end

disp('Case 3');
x0 = 1.45;
[x3, fval6, flag] = fsolve(@ej2, x0, options);

if flag == 1
    str = sprintf('El resultado es [%2.3f] partiendo del guess [x0=1.45]',x3);
    disp(str);
else
    str = sprintf('El resultado no converge partiendo del guess [x0=1.45]');
    disp(str);
end

disp('Case 4');
x0 = 3;
[x4, fval7, flag] = fsolve(@ej2, x0, options);

if flag == 1
    str=sprintf('El resultado es [%2.3f] partiendo del guess [x0=3]',x4);
    disp(str)
else
    str=sprintf('El resultado no converge partiendo del guess [x0=3]');
    disp(str)
end

% Grafica

X1 = linspace(1, 1.9999, 2500);
y1 = ((4*X1) - 7)./(X1 - 2);
X2 = linspace(2.0001, 3.5, 2500);
y2 = ((4*X2) - 7)./(X2 - 2);
X3 = linspace(1, 3.5, 5000);
T = zeros(5000, 1);
plot(X1, y1,'k-', ...
     X2, y2,'k-', ...
     X3, T, '--', ...
     1.625, 0, 'r*', ...
     1.875, 0, 'b*', ...
     1.45, 0, 'g*', ...
     3, 0, 'm*')

% Establece un rango para los ejes los primeros dos valores son
% para el eje "X" y los últimos dos son para el eje "Y".
axis([1, 3.6, -25, 25])
title('Ejercicio 1 inciso 2.');
xlabel('X');
ylabel('Y');

%
% Ejercicio 2
%

% Inciso 1

% Este comando sirve para que no se muetre en el command window
% "Equation solved..."
options = optimset('Display', 'off');
disp('Ejercicio 3, inciso 1.');
x0 = [1 1]';
[x, fval] = fsolve(@(x) f3i(x), x0, options);

% Este comando guarda en una variable el texto, el corchete
% dentro del texto [%2.3f,%2.3f] da la instrucción de sustituir
% las variables que aparecen al final del comando, en este caso
% x(1) y x(2).
str = sprintf('El resultado es [%2.3f, %2.3f] partiendo del guess [1,1]', ...
              x(1), x(2) );
disp(str);
