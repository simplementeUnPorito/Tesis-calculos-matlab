function p = params_martillo(varargin)
%PARAMS_MARTILLO Parametros del modelo leva-martillo (CAD de Fusion 360).
%
% Toda la geometria se mide del CAD con extraer_geometria_cad. Aca solo van
% las decisiones fisicas: material, velocidad de la leva, rigideces de
% contacto y el factor de escala CAD -> maquina real.
%
% Uso:
%   p = params_martillo();
%   p = params_martillo('rpm', 60, 'material', 'acero');
%   p = params_martillo('escala', 10);   % lleva el CAD 1:10 al diseno de campo
%
% Ver tambien: extraer_geometria_cad, construir_martillo_leva, simular_martillo

    q = inputParser();
    q.addParameter('rpm',      40,     @(x) isnumeric(x) && isscalar(x) && x > 0);
    q.addParameter('material', 'PLA',  @(s) ischar(s) || isstring(s));
    q.addParameter('escala',   1,      @(x) isnumeric(x) && isscalar(x) && x > 0);
    q.addParameter('t_fin',    [],     @(x) isempty(x) || isscalar(x));
    q.addParameter('recorte_nariz', Inf, @(x) isnumeric(x) && isscalar(x));
    q.addParameter('Geometria', [],    @(x) isempty(x) || isstruct(x));
    q.parse(varargin{:});
    o = q.Results;

    % --- geometria del CAD ---------------------------------------------------
    aqui = fileparts(mfilename('fullpath'));
    if ~isempty(o.Geometria)
        p.G = o.Geometria;
    else
        matg = fullfile(aqui, 'geometria_cad.mat');
        if isfile(matg)
            s = load(matg, 'G');  p.G = s.G;
        else
            p.G = extraer_geometria_cad('Guardar', true);
        end
    end
    p.assets = p.G.assets;

    % --- escala CAD -> maquina real -----------------------------------------
    % El CAD se dibujo para imprimir en 3D: las formas son las definitivas
    % pero el tamano es de maqueta. escala = 10 recupera el diseno de campo
    % (C = (250,400) mm, que es exactamente 10x el (25,40) del CAD).
    p.escala = o.escala;

    % --- material ------------------------------------------------------------
    p.material = char(o.material);
    switch lower(p.material)
        case 'pla',      p.rho = 1240;   % kg/m^3, impresion 3D
        case 'petg',     p.rho = 1270;
        case 'abs',      p.rho = 1040;
        case 'aluminio', p.rho = 2700;
        case 'acero',    p.rho = 7850;
        otherwise
            error('params_martillo:Material', 'Material desconocido: %s', p.material);
    end
    % la densidad no cambia con la escala, pero la masa va con escala^3
    p.rho_leva     = p.rho;
    p.rho_balancin = p.rho;
    p.rho_soporte  = p.rho;

    % --- accionamiento de la leva -------------------------------------------
    p.rpm   = o.rpm;
    p.w_cam = o.rpm * 2*pi/60;            % rad/s
    p.T_cam = 2*pi / p.w_cam;             % periodo de un golpe [s]

    % --- angulos --------------------------------------------------------------
    % theta = angulo del balancin. theta = 0 es EL IMPACTO: en esa pose la base
    % del mazo esta en z = -15 mm, el mismo nivel que la base de los soportes
    % (el piso). theta > 0 levanta el mazo.
    p.theta_impacto_deg = 0;
    p.theta_asm_deg     = p.G.theta_asm_deg;   % pose dibujada en el ensamblaje
    % Se arranca en la pose del ensamblaje: es la que dibujo el disenador, con
    % la leva ya apoyada en la cara del seguidor, asi que no hay un transitorio
    % inicial de penetracion.
    p.theta_ini_deg     = p.theta_asm_deg;

    % --- contacto leva-seguidor ---------------------------------------------
    % Penalizacion unilateral. La rigidez se elige para que la penetracion en
    % regimen sea del orden de 10 um con la carga de contacto tipica.
    m_bal = p.G.balancin.V_mm3 * 1e-9 * p.rho_balancin * p.escala^3;
    F_tip = m_bal * 9.80665 * 20;              % ~20x el peso, holgado
    p.pen_objetivo = 10e-6 * p.escala;         % m
    p.k_contacto   = F_tip / p.pen_objetivo;   % N/m
    p.zeta_contacto = 0.35;                     % amortiguamiento relativo
    p.c_contacto   = 2 * p.zeta_contacto * sqrt(p.k_contacto * m_bal);

    % --- recorte de la nariz de la leva (rediseno) ---------------------------
    % La leva del CAD llega a r = 45.1 mm, pero la cara del seguidor solo puede
    % ser empujada hasta rho_max = 42.17 mm (en theta = 58 deg, donde el punto
    % de contacto ya degenero sobre el pivote). Con la nariz entera el
    % mecanismo pasa sobre-centro y se traba. Este parametro recorta el radio
    % para evaluar el rediseno; Inf = la leva tal cual esta en el CAD.
    %   recorte_nariz en mm de radio, medido desde C.
    if isfinite(o.recorte_nariz)
        p.r_leva_max = o.recorte_nariz * 1e-3 * p.escala;   % m
    else
        p.r_leva_max = Inf;
    end
    p.recorte_nariz_mm = o.recorte_nariz;

    % --- yunque / piso (tope inferior del balancin) -------------------------
    % Se modela como el limite inferior de la junta revoluta en theta = 0.
    p.tope_rigidez  = 1e3 * p.escala^2;   % N*m/deg
    p.tope_amort    = 2.0 * p.escala^2;   % N*m/(deg/s)
    p.tope_ancho    = 0.05;               % deg, region de transicion
    p.tope_sup_deg  = 75;                 % tope superior de seguridad

    % --- friccion en los cojinetes ------------------------------------------
    p.amort_pivote = 1e-6 * p.escala^4;   % N*m/(rad/s)
    p.amort_leva   = 1e-6 * p.escala^4;

    % --- simulacion -----------------------------------------------------------
    if isempty(o.t_fin)
        p.t_fin = 2.2 * p.T_cam;          % dos golpes completos
    else
        p.t_fin = o.t_fin;
    end
    p.dt_max   = 2e-4;
    p.tau_filtro = 5e-4;                  % s, filtro del input de movimiento
    p.tol_rel  = 1e-6;
    p.gravedad = [0 0 -9.80665];          % z es la vertical del CAD

    % --- masas resultantes (informativo) ------------------------------------
    p.m_balancin = m_bal;
    p.m_leva     = p.G.leva.V_mm3 * 1e-9 * p.rho_leva * p.escala^3;
    p.m_soporte  = p.G.sop_der.V_mm3 * 1e-9 * p.rho_soporte * p.escala^3;

    % inercia del balancin respecto al pivote O (para chequeos analiticos)
    I_cm = p.G.balancin.I_rho_mm5(1,1) * 1e-15 * p.rho_balancin;   % kg m^2
    r_cm = norm(p.G.balancin.C_mm(2:3)) * 1e-3;                    % m
    m0   = p.G.balancin.V_mm3 * 1e-9 * p.rho_balancin;
    p.I_O = (I_cm + m0 * r_cm^2) * p.escala^5;
    p.S_bal = m0 * r_cm * p.escala^4;    % momento estatico [kg m]
end
