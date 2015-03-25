%c
% fiter.m    Programa para resolver una versión sencilla del problema del
%            planificador social usando el método de iteración de la función
%            de valor. La ecuación de Bellman a resolver es:
%
%            v(k) =  max {u(k) + beta v(k')}
%                    k'
%
function [kt, yt, it, ct] = v(alpha, beta, delta, A, ...
                              maxit, crit, T, kss, k0, malla)
    %
    % Mall para K
    %
    Kmin = 0.5*kss;
    Kmax = 1.2*kss;

    % Partición de la maya espaciado según:(Kmax-Kmin)/(malla-1)
    K = [Kmin:(Kmax - Kmin)/(malla - 1):Kmax]';

    % Guess inicial tamano 500x1
    V0 = 0*K;
    e = ones(malla, 1);

    % Matriz M es de 500 x 500 (ver slides)
    M = log(max((A*K.^alpha + (1 - delta)*K)*e' - e*K', 1e-8));

    %
    % Iteracion de la funcion de valor
    %
    figure(1)
    iconv = 0; % Indicadora, icon = 1 -> convergencia
    it = 1; % Iteraciones

    while (iconv == 0 & it < maxit)

        % Matriz 500x500, al transponer estamos sacando
        % el máximo de la primera fila inicial(pues
        % ya se convirtio en columna), la primera fila
        % indicaba F(K1, K1), F(K1, K2) ... F(K1, Kp).
        [V1, G] = max((M + beta*e*V0')');
        V1 = V1';
        G = G';

        if norm(V0 - V1) < crit
            iconv = 1;
        end

        %
        % Grafica
        %

        plot(K, V1);
        title(['Funcion de Valor en Iteración ', int2str(it), ...
               ' con Norma ', num2str(norm(V0 - V1))]);
        pause(0.01);
        V0 = V1; % Actualizacion.
        it = it + 1;
    end

    figure(2)

    % Dibujo de la función de politica del capital y la funcion identidad
    % (K,K). K(G) me indica dado el capital hoy cuál es el capital de mañana
    % y el capital de mañana se obtiene usando el vector de posicin óptimo
    % que maximiza la función valor.
    plot(K, K(G), K, K,':');
    title('Regla de Decisión Óptima');
    text(kss, kss, 'o  kss');

    % Simulación de las trayectorias óptimas:
    % - Capital (kt)
    % - Consumo (ct)
    % - Inversión (it)
    % - Producto (yt)
    kt = zeros(T + 1, 1);

    % Vector indicador de posiciones óptimas del capital que maximiza la
    % funcion valor en t+1
    ind = zeros(T, 1);

    % Recordar k0 = (2/3)*kss = 273.3665, esto nos dice cuando el vector
    % K es menor al k0 y su posición,en este caso para la posición 120 del
    % vector K tenemos un capital de 273.47< 273.365 esto es falso por tanto
    % coloca cero.
    [aux, ind(1)] = min(K < k0);

    % Vamos a partir de la posición K(120)= 273.47 una proxy del valor
    % inicial(k0) de (2/3)k de estado estacionario.

    kt(1)  = K(ind(1));
    for t = 1:T
        % G(ind(1)=120)=136 es la posicion óptima del capital que maximiza la
        % función valor que se usará para el siguiente periodo.
        ind(t + 1) = G(ind(t));

        % el capital de manana(k_2) óptimo va ser la posición 136 del vector
        % de capital K pues este maximiza la función valor el loop hace que
        % cada vez a medida que pasa el tiempo T el capital se aproxime a su
        % valor de estado estacionario.
        kt(t + 1)  = K(ind(t + 1));
    end

    yt = A*kt(1:T).^alpha; % Producto
    it = kt(2:T+1) - (1 - delta)*kt(1:T); % Inversion
    ct = yt - it; % Consumo