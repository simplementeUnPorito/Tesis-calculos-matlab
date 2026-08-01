function [Fy, Fz, T_bal, T_cam, delta, gap_brazo, py, pz] = ...
        contacto_leva_seguidor(theta, psi, w_bal, w_cam, k, c, escala, r_max)
%CONTACTO_LEVA_SEGUIDOR Contacto unilateral leva-seguidor, por penalizacion.
%
% El perfil de la leva NO es convexo (se desvia ~3.7 mm de su envolvente
% convexa justo en el escalon que produce la suelta), asi que el bloque
% Spatial Contact Force de Simscape -convexifica la geometria- no sirve:
% borraria el escalon. Aca el contacto se calcula analiticamente.
%
% La leva es ESTRELLADA respecto de su eje C, de modo que su contorno se
% describe con un radio univaluado r(alpha) y la deteccion cuesta O(1) por
% punto: se transforma cada punto del seguidor al marco de la leva y se compara
% su radio contra r(alpha).
%
% Entradas (SI, angulos en rad):
%   theta   angulo del balancin (0 = mazo en el piso, >0 = mazo arriba)
%   psi     angulo de la leva
%   w_bal   velocidad angular del balancin
%   w_cam   velocidad angular de la leva
%   k, c    rigidez [N/m] y amortiguamiento [N*s/m] del contacto
%   escala  factor CAD -> maquina real
%   r_max   recorte del radio de la nariz [m]. Inf = la leva del CAD tal cual.
%           Sirve para evaluar el rediseno: la leva del CAD llega a 45.1 mm
%           cuando la cara del seguidor solo puede alcanzar 42.17 mm, asi que
%           el mecanismo pasa sobre-centro y se traba en lugar de soltar.
%
% Salidas:
%   Fy, Fz     fuerza sobre el BALANCIN, en coordenadas del mundo [N]
%   T_bal      par de esa fuerza respecto al pivote O, alrededor de +X [N*m]
%   T_cam      par de la reaccion sobre la leva respecto a C, alrededor de +X
%   delta      penetracion maxima [m]  (<0 = separacion, no hay contacto)
%   gap_brazo  separacion minima leva-brazo del mazo [m] (chequeo de roce)
%   py, pz     punto de contacto en el mundo [m]  (diagnostico)

%#codegen

    d = coder.load('contacto_leva_datos.mat');
    alpha_t = d.alpha;                       % (NA,1) rad, uniforme en [0,2pi)
    r_t     = d.r_alpha * 1e-3 * escala;     % (NA,1) m
    if isfinite(r_max)
        r_t = min(r_t, r_max);               % nariz recortada (rediseno)
    end
    C       = d.C_cam(:)' * 1e-3 * escala;   % (1,2) eje de la leva en (y,z)
    Pc      = d.contacto * 1e-3 * escala;    % (n,2) region de contacto seguidor
    Pl      = d.poly_leva * 1e-3 * escala;   % (m,2) contorno de la leva

    NA = numel(alpha_t);
    dA = 2*pi / NA;

    ct = cos(theta);  st = sin(theta);
    cp = cos(psi);    sp = sin(psi);
    Rth = [ct -st; st ct];                   % marco balancin -> mundo
    Rps = [cp -sp; sp cp];                   % marco leva     -> mundo

    % ---- deteccion y fuerza: penalizacion DISTRIBUIDA ---------------------
    % Se suma la contribucion de TODOS los puntos penetrados, no solo la del
    % mas profundo. Tomar solo el maximo hace que la fuerza salte cada vez que
    % cambia cual es el punto mas profundo, y esa discontinuidad clava al
    % solver implicito de paso variable. Repartiendo la rigidez entre los N
    % puntos muestreados, cada uno aporta k/N*delta_i, que tiende a cero de
    % forma continua cuando el punto entra o sale del contacto; con
    % penetracion uniforme la fuerza total sigue siendo k*delta.
    N  = size(Pc, 1);
    ki = k / N;
    ci = c / N;

    Fy = 0;  Fz = 0;  T_bal = 0;  T_cam = 0;
    delta = -inf;                            % la maxima, solo como diagnostico
    py = 0;  pz = 0;
    sum_w = 0;

    for i = 1:N
        p = (Rth * Pc(i,:)')';               % mundo
        v = (Rps' * (p - C)')';              % marco de la leva
        rp = hypot(v(1), v(2));
        a  = atan2(v(2), v(1));
        if a < 0, a = a + 2*pi; end
        rc  = interp_periodica(r_t, a, dA, NA);
        del = rc - rp;
        if del > delta, delta = del; end

        if del <= 0
            continue
        end

        % normal exterior de la leva en a: para r(alpha),  n ~ r*er - r'*ea
        drd = deriv_periodica(r_t, a, dA, NA);
        er  = [cos(a); sin(a)];
        ea  = [-sin(a); cos(a)];
        nl  = rc*er - drd*ea;
        nl  = nl / max(norm(nl), eps);
        n   = Rps * nl;                      % normal en el mundo

        % velocidad relativa del punto de contacto a lo largo de la normal
        pv  = [p(1); p(2)];
        v_f = w_bal * [-pv(2); pv(1)];               % punto del seguidor
        rc2 = pv - C';
        v_c = w_cam * [-rc2(2); rc2(1)];             % punto de la leva
        vn  = (v_f - v_c)' * n;                      % >0 = se separan

        Fn = ki*del - ci*vn;
        if Fn < 0, Fn = 0; end                       % unilateral: solo empuja

        F  = Fn * n;
        Fy = Fy + F(1);
        Fz = Fz + F(2);
        T_bal = T_bal + pv(1)*F(2) - pv(2)*F(1);         % respecto a O, sobre +X
        T_cam = T_cam + rc2(1)*(-F(2)) - rc2(2)*(-F(1)); % reaccion, respecto a C

        % punto de contacto informado: promedio pesado por la penetracion
        py = py + del*p(1);  pz = pz + del*p(2);  sum_w = sum_w + del;
    end

    if sum_w > 0
        py = py / sum_w;  pz = pz / sum_w;
    end

    % ---- chequeo de roce contra el brazo del mazo --------------------------
    % puntos del contorno de la leva contra la recta del brazo del mazo
    % (el brazo va del pivote hacia +y en el marco del balancin)
    gap_brazo = inf;
    u = Rth * [1; 0];                        % direccion del brazo en el mundo
    n_arm = [-u(2); u(1)];                   % normal al brazo
    for j = 1:size(Pl, 1)
        q = (Rps * (Pl(j,:)' - C')) + C';    % mundo
        s = q' * u;                          % coordenada a lo largo del brazo
        if s > 0.01*escala && s < 0.090*escala
            gap_brazo = min(gap_brazo, abs(q' * n_arm) - 0.005*escala);
        end
    end

end

% ------------------------------------------------------------------------
function v = interp_periodica(tab, a, dA, NA)
%#codegen
    x  = a / dA;
    i0 = floor(x);
    f  = x - i0;
    i1 = mod(i0, NA) + 1;
    i2 = mod(i0 + 1, NA) + 1;
    v  = (1 - f) * tab(i1) + f * tab(i2);
end

function dv = deriv_periodica(tab, a, dA, NA)
%#codegen
    x  = a / dA;
    i0 = floor(x);
    im = mod(i0 - 1, NA) + 1;
    ip = mod(i0 + 1, NA) + 1;
    dv = (tab(ip) - tab(im)) / (2 * dA);
end
