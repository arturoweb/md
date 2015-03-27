%
% Programa para resolver una version sencilla del problema
% del planificador social usando el metodo Gauss-Seidel.
%
function [kt, yt, it, ct] = gs(alpha, ...
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

    % Partición de la malla
    K = [Kmin:(Kmax - Kmin)/(malla - 1):Kmax]';

    % Punto inicial
    V0 = 0*K;
    V1 = V0;
    n = length(K);
    e = ones(malla, 1);

    % Matriz M es de 500 x 500
    M = log(max((A*K.^alpha + (1 - delta)*K)*e' - e*K', 1e-8));

    figure(1)
    iconv = 0; % Indicadora, iconv = 1 -> convergencia
    it = 1;    % Iteraciones

    while (iconv == 0 & it < maxit)

        %
        % Gaus-Seidel
        %
        for i = 1 : n
            V1(i) = V0(i);
            V1(i) = V1(i) - M(i, 1:i - 1) * V1(1:i - 1);
            V1(i) = V1(i) - M(i, i + 1:n) * V0(i + 1:n);
            V1(i) = V1(i)/M(i, i);
        end

        if norm(V0 - V1) < crit
            iconv = 1;
        end

        V0 = V1; % Actualizacion
        it = it + 1;
    end

    % Simulacion de trayectorias optimas
    kt = zeros(T + 1, 1);
    ind = zeros(T, 1);

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