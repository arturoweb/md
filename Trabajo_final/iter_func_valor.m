%
% Método de iteración de la función valor
%
% Trabajo final - Macroeconomía Dinámica
% Profesor: Carlos Urrutia
% ITAM, 2015
%
% Equipo:
% Omar Trejo, 119711
% Alejandro Cerecero, 000000
% Arturo Reynoso, 000000
%
function [kt, yt, it, ct] = iter_func_valor(alpha, ...
                                            beta, ...
                                            delta, ...
                                            A, ...
                                            max_iter, ...
                                            crit, ...
                                            T, ...
                                            kss, ...
                                            k0, ...
                                            malla)
    % TODO: (otn) ¿Por qué estos valores?
    
    % Malla para K
    k_min = 0.8*kss;
    k_max = 1.5*kss;

    % Partición de la malla
    K = [k_min:(k_max - k_min)/(malla - 1):k_max]';

    % Punto inicial
    V0 = 0*K;
    e = ones(malla, 1);

    % Matriz M
    M = log(max((A*K.^alpha + (1 - delta)*K)*e' - e*K', 1e-8));

    %
    % Iteracion de la funcion de valor
    %
    figure(1)
    terminado = 0;
    iter = 1;

    while terminado == 0 & iter < max_iter

        % Al transponer estamos sacando
        % el máximo de la primera fila inicial(pues
        % ya se convirtio en columna), la primera fila
        % indicaba F(K1, K1), F(K1, K2) ... F(K1, Kp).
        [V1, G] = max((M + beta*e*V0')');
        V1 = V1';
        G = G';

        if norm(V0 - V1) < crit
            terminado = 1;
        end

        % Gráfica
        plot(K, V1);
        title(['Valor en iteracion ', int2str(iter), ...
               ' con norma ', num2str(norm(V0 - V1))]);
        pause(0.01);
        V0 = V1;
        iter = iter + 1;
    end

    figure(2)
    % Gráfica de la función de politica de
    % capital y la funcion identidad (K,K)
    plot(K, K(G), K, K,':');
    title('Regla de Decisión Óptima');
    text(kss, kss, 'o  kss');

    % Simulacion de trayectorias optimas
    kt = zeros(T + 1, 1);
    ind = zeros(T, 1);

    % Recordar k0 = (2/3)*kss = 273.3665, esto nos dice cuando 
    % el vector K es menor al k0 y su posicion, en este caso 
    % para la posición 120 del vector K tenemos un capital de 
    % 273.47 < 273.365 esto es falso por tanto coloca cero.
    
    [aux, ind(1)] = min(K < k0);
    kt(1) = K(ind(1));
    
    for t = 1:T
        ind(t + 1) = G(ind(t));
        kt(t + 1)  = K(ind(t + 1));
    end

    yt = A*kt(1:T).^alpha;                % Producto
    it = kt(2:T+1) - (1 - delta)*kt(1:T); % Inversion
    ct = yt - it;                         % Consumo
end
