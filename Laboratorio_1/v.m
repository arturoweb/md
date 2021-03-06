%
% Programa para resolver una version sencilla del problema del
% planificador social usando el metodo de iteracion de la funcion
% de valor. La ecuacion de Bellman a resolver es:
%
% v(k) =  max {u(k) + beta v(k')}
% k'
%
function [kt, yt, it, ct] = v(alpha, ...
                              beta, ...
                              delta, ...
                              A, ...
                              maxit, ...
                              crit, ...
                              T, ...
                              kss, ...
                              k0, ...
                              malla)
    % Malla para K
    Kmin = 0.5*kss;
    Kmax = 1.2*kss;

    % Partici�n de la malla
    K = [Kmin:(Kmax - Kmin)/(malla - 1):Kmax]';

    % Punto inicial
    V0 = 0*K;
    e = ones(malla, 1);

    % Matriz M es de 500 x 500
    M = log(max((A*K.^alpha + (1 - delta)*K)*e' - e*K', 1e-8));

    %
    % Iteracion de la funcion de valor
    %
    figure(1)
    iconv = 0; % Indicadora, iconv = 1 -> convergencia
    it = 1;    % Iteraciones

    while (iconv == 0 & it < maxit)

        % Matriz 500x500, al transponer estamos sacando
        % el m�ximo de la primera fila inicial(pues
        % ya se convirtio en columna), la primera fila
        % indicaba F(K1, K1), F(K1, K2) ... F(K1, Kp).
        [V1, G] = max((M + beta*e*V0')');
        V1 = V1';
        G = G';

        if norm(V0 - V1) < crit
            iconv = 1;
        end

        % Grafica
        plot(K, V1);
        title(['Funcion de Valor en Iteraci�n ', int2str(it), ...
               ' con Norma ', num2str(norm(V0 - V1))]);
        pause(0.01);
        V0 = V1;
        it = it + 1;
    end

    figure(2)

    % Dibujo de la funci�n de politica del capital y la funcion identidad
    % (K,K). K(G) me indica dado el capital hoy cual es el capital de manana
    % y el capital de manana se obtiene usando el vector de posicion optimo
    % que maximiza la funcion valor.
    plot(K, K(G), K, K,':');
    title('Regla de Decisi�n �ptima');
    text(kss, kss, 'o  kss');

    % Simulacion de trayectorias optimas
    kt = zeros(T + 1, 1);
    ind = zeros(T, 1);

    % Recordar k0 = (2/3)*kss = 273.3665, esto nos dice cuando el vector
    % K es menor al k0 y su posicion, en este caso para la posici�n 120 del
    % vector K tenemos un capital de 273.47 < 273.365 esto es falso por tanto
    % coloca cero.
    [aux, ind(1)] = min(K < k0);

    % Vamos a partir de la posici�n K(120)= 273.47 una proxy del valor
    % inicial(k0) de (2/3)k de estado estacionario.
    kt(1) = K(ind(1));
    for t = 1:T
        % G(ind(1)=120)=136 es la posicion optima del capital que maximiza la
        % funcion valor que se usara para el siguiente periodo.
        ind(t + 1) = G(ind(t));

        % El capital de manana(k_2) optimo va ser la posicion 136 del vector
        % de capital K pues este maximiza la funcion valor el loop hace que
        % cada vez a medida que pasa el tiempo T el capital se aproxime a su
        % valor de estado estacionario.
        kt(t + 1)  = K(ind(t + 1));
    end

    yt = A*kt(1:T).^alpha;                % Producto
    it = kt(2:T+1) - (1 - delta)*kt(1:T); % Inversion
    ct = yt - it;                         % Consumo
end
