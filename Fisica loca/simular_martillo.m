function R = simular_martillo(varargin)
%SIMULAR_MARTILLO Corre el modelo del martillo de leva y reporta los resultados.
%
% Responde las cuatro preguntas de diseno:
%   1. impacto      velocidad y energia del mazo al golpear, y repetibilidad
%   2. suelta       si la leva realmente suelta y si roza durante la caida
%   3. motor        par requerido en la leva y limites de rpm
%   4. esfuerzos    reacciones en el pivote O y en el eje de la leva C
%
% Las reacciones y el par del motor se calculan en post-proceso con la dinamica
% del cuerpo rigido (el mecanismo es plano y de 1 grado de libertad, asi que es
% exacto) en vez de con puertos de sensado de las juntas.
%
% Uso:
%   R = simular_martillo();                       % la leva tal cual el CAD
%   R = simular_martillo('recorte_nariz', 41);    % nariz recortada a r=41 mm
%   R = simular_martillo('rpm', 60, 'Figuras', true);
%
% Ver tambien: params_martillo, construir_martillo_leva

    q = inputParser();
    q.KeepUnmatched = true;
    q.addParameter('Figuras', true, @islogical);
    q.addParameter('Reconstruir', false, @islogical);
    q.parse(varargin{:});
    opt = q.Results;
    resto = q.Unmatched;

    aqui = fileparts(mfilename('fullpath'));
    old = cd(aqui);  limpiar = onCleanup(@() cd(old));

    args = reshape([fieldnames(resto)'; struct2cell(resto)'], 1, []);
    p = params_martillo(args{:});

    mdl = 'martillo_leva';
    if opt.Reconstruir || ~isfile(fullfile(aqui, [mdl '.slx']))
        construir_martillo_leva(p);
    end
    if ~bdIsLoaded(mdl), load_system(mdl); end
    set_param(mdl, 'SimscapeLogLimitData', 'off');

    assignin('base', 'p', p);
    fprintf('corriendo %s: %.1f rpm, %s, escala %g, t_fin %.2f s ...\n', ...
            mdl, p.rpm, p.material, p.escala, p.t_fin);
    out = sim(mdl, 'ReturnWorkspaceOutputs', 'on');

    R = analizar(out, p);
    if opt.Figuras
        figuras(R, p, aqui);
    end
    reporte(R, p);
end

% =========================================================================
function R = analizar(out, p)
    sl = out.simlog;
    G  = p.G;

    R.t     = sl.Pivote.Rz.q.series.time;
    R.theta = sl.Pivote.Rz.q.series.values('rad');
    R.omega = sl.Pivote.Rz.w.series.values('rad/s');
    R.psi   = interp1(sl.Leva_j.Rz.q.series.time, ...
                      sl.Leva_j.Rz.q.series.values('rad'), R.t, 'linear', 'extrap');

    % diagnostico del contacto, remuestreado sobre la malla de simlog
    D = out.diag.signals.values;  td = out.diag.time;
    R.delta = interp1(td, D(:,1), R.t, 'linear', 'extrap');
    R.gap   = interp1(td, D(:,2), R.t, 'linear', 'extrap');
    R.pc    = [interp1(td, D(:,3), R.t, 'linear', 'extrap'), ...
               interp1(td, D(:,4), R.t, 'linear', 'extrap')];

    % aceleracion angular por diferencias centradas
    R.alpha = gradiente(R.omega, R.t);
    R.contacto = R.delta > 0;

    % --- fuerza de contacto reconstruida ---------------------------------
    % Se recalcula con la misma funcion del modelo, para tener Fy, Fz y los
    % pares sin depender de mas puertos de sensado.
    n = numel(R.t);
    R.F = zeros(n, 2);  R.T_bal = zeros(n,1);  R.T_cam = zeros(n,1);
    for i = 1:n
        [fy, fz, tb, tc] = contacto_leva_seguidor(R.theta(i), R.psi(i), ...
            R.omega(i), p.w_cam, p.k_contacto, p.c_contacto, p.escala, p.r_leva_max);
        R.F(i,:) = [fy fz];  R.T_bal(i) = tb;  R.T_cam(i) = tc;
    end

    % --- geometria de masas -----------------------------------------------
    esc  = p.escala;
    m_b  = p.m_balancin;
    c_b  = G.balancin.C_mm(2:3)' * 1e-3 * esc;      % centroide en el marco propio
    m_l  = p.m_leva;
    c_l  = G.leva.C_mm(2:3)' * 1e-3 * esc;
    C_c  = G.C' * 1e-3 * esc;
    g    = [0; p.gravedad(3)];

    % --- reaccion en el pivote O -------------------------------------------
    % m*a_cm = F_contacto + m*g + R_O
    R.R_O = zeros(n, 2);
    for i = 1:n
        th = R.theta(i);
        Rot = [cos(th) -sin(th); sin(th) cos(th)];
        r   = Rot * c_b;                              % centroide en el mundo
        a   = R.alpha(i) * [-r(2); r(1)] - R.omega(i)^2 * r;
        R.R_O(i,:) = (m_b*a - R.F(i,:)' - m_b*g)';
    end
    R.R_O_mod = vecnorm(R.R_O, 2, 2);

    % --- reaccion en el eje de la leva C ----------------------------------
    % la leva gira a velocidad constante: su centroide describe un circulo
    R.R_C = zeros(n, 2);
    for i = 1:n
        ps = R.psi(i);
        Rot = [cos(ps) -sin(ps); sin(ps) cos(ps)];
        r   = C_c + Rot * (c_l - C_c);
        a   = -p.w_cam^2 * (r - C_c);                 % aceleracion centripeta
        R.R_C(i,:) = (m_l*a + R.F(i,:)' - m_l*g)';    % -F es la reaccion
    end
    R.R_C_mod = vecnorm(R.R_C, 2, 2);

    % --- par del motor -----------------------------------------------------
    % leva a velocidad constante -> el par equilibra el contacto y la gravedad
    R.T_motor = zeros(n, 1);
    for i = 1:n
        ps = R.psi(i);
        Rot = [cos(ps) -sin(ps); sin(ps) cos(ps)];
        r   = Rot * (c_l - C_c);                      % C -> centroide
        T_grav = r(1)*m_l*g(2) - r(2)*m_l*g(1);
        R.T_motor(i) = -R.T_cam(i) - T_grav;
    end

    % --- eventos: suelta e impacto ----------------------------------------
    R.I = eventos(R, p);
end

% =========================================================================
function I = eventos(R, p)
    th = R.theta;  w = R.omega;  t = R.t;
    tope = deg2rad(p.tope_sup_deg);

    % impactos: cruce descendente de theta por 0.2 deg con velocidad negativa
    u = deg2rad(0.2);
    ii = find(th(1:end-1) > u & th(2:end) <= u & w(2:end) < 0);
    I.t_impacto = t(ii);
    I.w_impacto = w(ii);
    Lh = p.G.balancin.y_mazo * 1e-3 * p.escala;    % brazo del mazo
    I.Lh = Lh;
    I.v_impacto = abs(w(ii)) * Lh;
    I.E_impacto = 0.5 * p.I_O * w(ii).^2;

    % alzada maxima de cada ciclo y si toco el tope de seguridad
    I.theta_max = max(th);
    I.toca_tope = I.theta_max >= tope - deg2rad(0.5);

    % suelta: ultimo instante con contacto antes de cada impacto
    I.t_suelta = nan(size(ii));
    I.w_suelta = nan(size(ii));
    I.theta_suelta = nan(size(ii));
    for j = 1:numel(ii)
        k = find(R.contacto(1:ii(j)), 1, 'last');
        if ~isempty(k)
            I.t_suelta(j) = t(k);
            I.w_suelta(j) = w(k);
            I.theta_suelta(j) = th(k);
        end
    end
    I.t_caida = I.t_impacto - I.t_suelta;

    % interferencia durante la caida (fuera del contacto)
    libre = ~R.contacto;
    I.gap_min_libre = min(R.gap(libre));
    I.penetracion_max = max(R.delta);

    % El piso se modela como par de tope de la junta, no como una fuerza
    % externa sobre el mazo. Mientras el balancin esta apoyado contra el tope,
    % la reaccion calculada en O absorbe tambien la fuerza del piso y no
    % representa la carga del cojinete: esas muestras se marcan aparte.
    I.en_tope = th < deg2rad(0.05) | th > tope - deg2rad(0.05);
end

% =========================================================================
function d = gradiente(y, t)
    d = zeros(size(y));
    d(2:end-1) = (y(3:end) - y(1:end-2)) ./ (t(3:end) - t(1:end-2));
    if numel(y) > 1
        d(1) = (y(2)-y(1))/(t(2)-t(1));
        d(end) = (y(end)-y(end-1))/(t(end)-t(end-1));
    end
end

% =========================================================================
function figuras(R, p, aqui)
    I = R.I;
    f = figure('Position', [50 50 1350 900], 'Visible', 'off');
    tl = tiledlayout(f, 4, 2, 'TileSpacing','compact', 'Padding','compact');
    title(tl, sprintf('martillo de leva  -  %.0f rpm, %s, escala %g, nariz %s', ...
        p.rpm, p.material, p.escala, recorte_txt(p)), 'FontWeight','bold');

    nexttile([1 2]);
    plot(R.t, rad2deg(R.theta), 'LineWidth', 1.2); grid on; hold on
    yline(0, 'k--'); yline(p.tope_sup_deg, 'r--', 'tope de seguridad');
    plot(I.t_impacto, zeros(size(I.t_impacto)), 'rv', 'MarkerFaceColor','r');
    plot(I.t_suelta, rad2deg(I.theta_suelta), 'go', 'MarkerFaceColor','g');
    ylabel('\theta [deg]'); legend('\theta','','','impacto','suelta','Location','best');
    title('angulo del balancin (0 = mazo en el piso)');

    nexttile; plot(R.t, R.omega, 'LineWidth', 1.1); grid on
    ylabel('\omega [rad/s]'); title('velocidad angular');

    nexttile; plot(R.t, R.delta*1e3, 'LineWidth', 1.1); grid on; hold on
    yline(0,'k--'); ylabel('\delta [mm]'); title('penetracion del contacto');

    nexttile; plot(R.t, vecnorm(R.F,2,2), 'LineWidth', 1.1); grid on
    ylabel('|F| [N]'); title('fuerza de contacto leva-seguidor');

    nexttile; plot(R.t, R.T_motor*1e3, 'LineWidth', 1.1); grid on
    ylabel('T_{motor} [mN\cdotm]'); title('par en el eje de la leva');

    nexttile; plot(R.t, R.R_O_mod, 'LineWidth', 1.1); grid on
    ylabel('|R_O| [N]'); xlabel('t [s]'); title('reaccion en el pivote O');

    nexttile; plot(R.t, R.R_C_mod, 'LineWidth', 1.1); grid on
    ylabel('|R_C| [N]'); xlabel('t [s]'); title('reaccion en el eje de la leva');

    dest = fullfile(aqui, 'martillo_leva_curvas.png');
    exportgraphics(f, dest, 'Resolution', 120);
    close(f);
    fprintf('-> %s\n', dest);
end

function s = recorte_txt(p)
    if isfinite(p.recorte_nariz_mm)
        s = sprintf('recortada a %.1f mm', p.recorte_nariz_mm);
    else
        s = 'del CAD (sin recortar)';
    end
end

% =========================================================================
function reporte(R, p)
    I = R.I;  G = p.G;
    L = @(varargin) fprintf(varargin{:});
    L('\n%s\n', repmat('=', 1, 72));
    L('  MARTILLO DE LEVA  -  %.0f rpm, %s, escala %g, nariz %s\n', ...
      p.rpm, p.material, p.escala, recorte_txt(p));
    L('%s\n', repmat('=', 1, 72));

    L('\nMASAS\n');
    L('  balancin (brazo + mazo + seguidor) : %8.2f g\n', p.m_balancin*1e3);
    L('  leva                               : %8.2f g\n', p.m_leva*1e3);
    L('  I_O (balancin en el pivote)        : %10.4e kg m^2\n', p.I_O);

    L('\n1. IMPACTO\n');
    if isempty(I.t_impacto)
        L('  el mazo NO llega al piso en la ventana simulada\n');
    else
        L('  golpes detectados: %d\n', numel(I.t_impacto));
        L('    %-8s %-12s %-12s %-12s\n', 't [s]', 'w [rad/s]', 'v [m/s]', 'E [mJ]');
        for j = 1:numel(I.t_impacto)
            L('    %-8.4f %-12.3f %-12.4f %-12.3f\n', ...
              I.t_impacto(j), I.w_impacto(j), I.v_impacto(j), I.E_impacto(j)*1e3);
        end
        if numel(I.v_impacto) > 1
            L('  repetibilidad de v: media %.4f m/s, desvio %.5f m/s (%.2f%%)\n', ...
              mean(I.v_impacto), std(I.v_impacto), ...
              100*std(I.v_impacto)/mean(I.v_impacto));
        end
        L('  brazo del mazo usado: %.1f mm\n', I.Lh*1e3);
    end

    L('\n2. SUELTA Y NO-INTERFERENCIA\n');
    L('  alzada maxima          : %.2f deg', rad2deg(I.theta_max));
    if I.toca_tope
        L('   <-- LLEGA AL TOPE DE SEGURIDAD (%g deg)\n', p.tope_sup_deg);
    else
        L('\n');
    end
    if ~isempty(I.theta_suelta) && any(~isnan(I.theta_suelta))
        L('  suelta en theta        : %.2f deg\n', rad2deg(mean(I.theta_suelta,'omitnan')));
        L('  velocidad en la suelta : %.4f rad/s (deberia ser ~0)\n', ...
          mean(I.w_suelta,'omitnan'));
        L('  tiempo de caida        : %.4f s\n', mean(I.t_caida,'omitnan'));
    end
    L('  penetracion maxima     : %.3f mm\n', I.penetracion_max*1e3);
    L('  separacion minima leva-brazo fuera del contacto: %.3f mm\n', ...
      I.gap_min_libre*1e3);

    % veredicto sobre el sobre-centro
    rho_max = G.rho_max_alcanzable;
    r_leva  = G.leva.r_max;
    L('\n  chequeo cinematico de la leva:\n');
    L('    radio maximo de la leva          : %.2f mm\n', r_leva);
    L('    rho maximo alcanzable por la cara: %.2f mm (en theta = %.1f deg)\n', ...
      rho_max, rad2deg(G.theta_rho_max));
    if r_leva > rho_max
        L('    [FALLA] la nariz excede en %.2f mm: el punto de contacto migra\n', ...
          r_leva - rho_max);
        L('            hacia el pivote y el mecanismo pasa SOBRE-CENTRO.\n');
        L('            Recorte sugerido de la nariz: r <= %.1f mm\n', ...
          0.97*rho_max);
    else
        L('    [OK] la nariz entra en el rango util de la cara\n');
    end

    L('\n3. MOTOR\n');
    L('  par en la leva: medio %.3f mN m, pico %.3f mN m\n', ...
      mean(abs(R.T_motor))*1e3, max(abs(R.T_motor))*1e3);
    L('  potencia pico : %.4f W\n', max(abs(R.T_motor))*abs(p.w_cam));
    if ~isempty(I.t_caida) && any(~isnan(I.t_caida))
        tc = mean(I.t_caida,'omitnan');
        L('  tiempo de caida %.4f s -> la leva no debe volver antes:\n', tc);
        L('    rpm maxima admisible ~ %.1f rpm\n', 60/(tc/0.55));
    end

    L('\n4. ESFUERZOS EN EJES\n');
    fuera = ~I.en_tope;
    L('  fuera de los topes (carga real de cojinete):\n');
    L('    |R_O| pivote   : medio %.3f N, pico %.3f N\n', ...
      mean(R.R_O_mod(fuera)), max(R.R_O_mod(fuera)));
    L('    |R_C| eje leva : medio %.3f N, pico %.3f N\n', ...
      mean(R.R_C_mod(fuera)), max(R.R_C_mod(fuera)));
    L('  peso del balancin, como referencia: %.3f N\n', p.m_balancin*9.80665);
    L('  NOTA: el piso esta modelado como par de tope de la junta, no como una\n');
    L('  fuerza externa sobre el mazo. Durante el apoyo contra el tope la\n');
    L('  reaccion en O absorbe tambien la fuerza del piso, asi que esas\n');
    L('  muestras (%.1f%% del tiempo) se excluyen del maximo de arriba. Para\n', ...
      100*mean(I.en_tope));
    L('  dimensionar el eje en el instante del golpe hace falta modelar el\n');
    L('  yunque como contacto explicito, no como limite de junta.\n');
    L('\n%s\n', repmat('=', 1, 72));
end
