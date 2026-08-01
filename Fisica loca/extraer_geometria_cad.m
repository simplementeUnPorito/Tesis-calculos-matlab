function G = extraer_geometria_cad(varargin)
%EXTRAER_GEOMETRIA_CAD Mide el mecanismo leva-martillo desde el CAD de Fusion 360.
%
% Lee las mallas STL del ensamblaje y devuelve la geometria y las propiedades
% masicas que necesita el modelo de Simscape Multibody. Todo se mide, nada se
% asume: los ejes de pivote salen de los bujes cilindricos y el perfil de la
% leva de un corte de la malla.
%
% Sistema de coordenadas del CAD (mm):
%   - el mecanismo es plano y vive en el plano YZ
%   - z es la vertical (hacia arriba), x es el eje de los pernos
%   - O = (y,z) = (0,0)   pivote del balancin
%   - C = (y,z) = (25,40) eje de la leva
%
% Uso:
%   G = extraer_geometria_cad();                  % usa Assets/Fusion 360
%   G = extraer_geometria_cad('Guardar', true);   % guarda geometria_cad.mat
%
% Ver tambien: params_martillo, construir_martillo_leva

    p = inputParser();
    p.addParameter('Assets', '', @(s) ischar(s) || isstring(s));
    p.addParameter('Guardar', true, @islogical);
    p.parse(varargin{:});
    opt = p.Results;

    assets = char(opt.Assets);
    if isempty(assets)
        aqui = fileparts(mfilename('fullpath'));
        % .../src/calculos_modelados/matlab/Fisica loca -> raiz del superproyecto
        assets = fullfile(aqui, '..', '..', '..', '..', 'Assets', 'Fusion 360');
    end
    if ~isfolder(assets)
        error('extraer_geometria_cad:SinAssets', ...
              'No encuentro la carpeta de CAD: %s', assets);
    end
    G.assets = char(java.io.File(assets).getCanonicalPath());

    % --- ejes del mecanismo, medidos de los bujes (ver README) --------------
    G.O = [0 0];        % pivote del balancin, en (y,z) [mm]
    G.C = [25 40];      % eje de la leva,     en (y,z) [mm]
    G.eje = [1 0 0];    % los dos giran alrededor de +X

    % --- mallas -------------------------------------------------------------
    piezas = struct( ...
        'balancin', 'Mazo.stl', ...
        'leva',     'Leva.stl', ...
        'sop_der',  'Base derecha.stl', ...
        'sop_izq',  'Base izquierda.stl');
    campos = fieldnames(piezas);

    for k = 1:numel(campos)
        f = fullfile(G.assets, piezas.(campos{k}));
        if ~isfile(f)
            error('extraer_geometria_cad:SinSTL', 'Falta el STL: %s', f);
        end
        T = stlread(f);
        V = T.Points;                 % (n,3) mm
        F = T.ConnectivityList;
        [Vol, Cm, I] = props_masicas(V, F);
        G.(campos{k}) = struct( ...
            'archivo', piezas.(campos{k}), ...
            'V_mm3',   Vol, ...
            'C_mm',    Cm, ...          % centroide, coords del archivo
            'I_rho_mm5', I, ...         % inercia en el centroide, densidad 1
            'bb_min',  min(V, [], 1), ...
            'bb_max',  max(V, [], 1), ...
            'V_puntos', V, 'F_caras', F);
    end

    % --- perfil de la leva: corte del disco en x = 0 ------------------------
    % El contorno NO se encadena en un poligono: encadenar segmentos sueltos es
    % fragil y un bucle mal cerrado mete una cuerda recta que en la simulacion
    % aparece como un flanco empinado inexistente. Como la leva es estrellada
    % respecto de C, alcanza con intersecar rayos desde C contra los segmentos
    % crudos, y el poligono (solo para graficar) se obtiene ordenando por angulo.
    L = G.leva;
    segsL = corte_yz(L.V_puntos, L.F_caras, 0.0);
    G.leva.n_segmentos = size(segsL, 1);
    polyL = ordenar_angular(segsL, G.C);
    G.leva.poligono = polyL;                    % (n,2) en (y,z) [mm]
    rL = vecnorm(polyL - G.C, 2, 2);
    G.leva.r_min = min(rL);
    G.leva.r_max = max(rL);
    % espesor axial del disco (los triangulos lejos del eje)
    ctro = squeeze(mean(reshape(L.V_puntos(L.F_caras', :), 3, [], 3), 1));
    rad  = vecnorm(ctro(:, 2:3) - G.C, 2, 2);
    xd   = ctro(rad > 20, 1);
    G.leva.x_disco = [min(xd) max(xd)];

    % perfil radial r(alpha) respecto de C: la leva es estrellada en C, asi que
    % r(alpha) es univaluada y sirve como tabla de consulta para el contacto
    NA = 1440;
    G.leva.alpha = linspace(0, 2*pi, NA + 1)';  G.leva.alpha(end) = [];
    G.leva.r_alpha = radio_por_rayo(segsL, G.C, G.leva.alpha);

    % Validacion del perfil. Un contorno mal extraido aparece como un salto
    % radial grande, y en la simulacion eso golpea al seguidor como si la leva
    % tuviera un flanco vertical que no existe. Se corta aca antes de simular.
    if any(~isfinite(G.leva.r_alpha))
        error('extraer_geometria_cad:PerfilConHuecos', ...
            ['r(alpha) tiene %d direcciones sin interseccion: el corte de la ' ...
             'leva no cerro o el perfil no es estrellado respecto de C.'], ...
            sum(~isfinite(G.leva.r_alpha)));
    end
    dr_max = max(abs(diff([G.leva.r_alpha; G.leva.r_alpha(1)])));
    paso_deg = 360 / NA;
    G.leva.pendiente_max = dr_max / paso_deg;      % mm de radio por grado
    if G.leva.pendiente_max > 3.0
        warning('extraer_geometria_cad:PerfilAbrupto', ...
            ['r(alpha) salta hasta %.2f mm/deg. Si el perfil real es suave, ' ...
             'la extraccion metio una cuerda falsa y la simulacion va a dar ' ...
             'un golpe inexistente.'], G.leva.pendiente_max);
    end

    % --- balancin: corte en x = 0, cara del seguidor y punta ---------------
    % Aca solo se necesita el conjunto de puntos del contorno (para medir la
    % cara, la punta y muestrear la region de contacto), no un poligono
    % ordenado: el balancin es una L y no es estrellado respecto del pivote.
    B = G.balancin;
    segsB = corte_yz(B.V_puntos, B.F_caras, 0.0);
    polyB = unique(round(reshape(permute(segsB, [1 2 3]), [], 2), 6), 'rows');
    G.balancin.poligono = polyB;
    G.balancin.z_punta  = max(polyB(:, 2));     % punta del brazo seguidor
    cara = polyB(polyB(:, 2) > 10 & abs(polyB(:, 1)) < 12, :);
    G.balancin.d_cara   = max(cara(:, 1));      % cara plana del seguidor en y=+d
    G.balancin.y_mazo   = max(polyB(:, 1));
    G.balancin.z_mazo   = min(polyB(:, 2));

    % region de contacto del seguidor: cara plana + punta redondeada
    sel = polyB(:, 2) > 2.0 & polyB(:, 1) > -1.0;
    cont = polyB(sel, :);
    [~, ord] = sort(cont(:, 2));
    G.balancin.contacto = cont(ord, :);

    % --- cinematica derivada ------------------------------------------------
    % distancia del eje de la leva a la cara del seguidor cuando el balancin
    % esta en theta:   rho(theta) = Cy*cos(theta) + Cz*sin(theta) - d_cara
    % (se deja como coeficientes, no como function handle: asi el .mat queda
    %  con puros datos numericos y lo puede leer coder.load)
    G.rho_coef = [G.C(1) G.C(2) -G.balancin.d_cara];
    G.rho_max_alcanzable = hypot(G.C(1), G.C(2)) - G.balancin.d_cara;
    G.theta_rho_max = atan2(G.C(2), G.C(1));
    G.theta_asm_deg = 15.897;    % pose del balancin en el ensamblaje STEP

    % --- limpieza: las mallas crudas no van al .mat -------------------------
    for k = 1:numel(campos)
        G.(campos{k}) = rmfield(G.(campos{k}), {'V_puntos', 'F_caras'});
    end

    informe(G);

    if opt.Guardar
        aqui = fileparts(mfilename('fullpath'));

        % (1) geometria completa, para los scripts de analisis
        d1 = fullfile(aqui, 'geometria_cad.mat');
        save(d1, 'G');

        % (2) solo los arreglos numericos del contacto. Va aparte porque
        % coder.load -que usa el bloque MATLAB Function- carga TODAS las
        % variables del archivo y no tolera structs anidados grandes.
        alpha     = G.leva.alpha;            %#ok<NASGU>
        r_alpha   = G.leva.r_alpha;          %#ok<NASGU>
        C_cam     = G.C;                     %#ok<NASGU>
        contacto  = G.balancin.contacto;     %#ok<NASGU>
        poly_leva = G.leva.poligono;         %#ok<NASGU>
        d2 = fullfile(aqui, 'contacto_leva_datos.mat');
        save(d2, 'alpha', 'r_alpha', 'C_cam', 'contacto', 'poly_leva');

        fprintf('\n-> %s\n-> %s\n', d1, d2);
    end
end

% =========================================================================
function [V, Cm, I] = props_masicas(P, F)
%PROPS_MASICAS Volumen, centroide e inercia de una malla triangular cerrada.
% Densidad unitaria; V en mm^3, Cm en mm, I en mm^5 (multiplicar por rho).
    a = P(F(:,1), :);  b = P(F(:,2), :);  c = P(F(:,3), :);
    v6 = sum(a .* cross(b, c, 2), 2);          % 6*volumen con signo del tet
    V  = sum(v6) / 6;
    Cm = sum(((a + b + c) / 4) .* (v6 / 6), 1) / V;

    % segundo momento respecto al origen, por tetraedro (origen,a,b,c)
    S = zeros(3);
    w = v6 / 6;
    for i = 1:size(a, 1)
        Pt = [0 0 0; a(i,:); b(i,:); c(i,:)];
        s  = zeros(3);
        for m = 1:4
            for n = 1:4
                s = s + Pt(m,:)' * Pt(n,:) * (1 + (m == n));
            end
        end
        S = S + s * (w(i) / 20);
    end
    I = trace(S) * eye(3) - S;
    I = I - V * (dot(Cm, Cm) * eye(3) - Cm' * Cm);   % trasladar al centroide
end

% =========================================================================
function segs = corte_yz(P, F, x0)
%CORTE_YZ Contorno (y,z) de la interseccion de la malla con el plano x = x0.
    x = P(:, 1);
    xf = x(F);                                  % (nf,3)
    cruza = ~(all(xf > x0, 2) | all(xf < x0, 2));
    Fc = F(cruza, :);
    segs = zeros(0, 2, 2);
    for i = 1:size(Fc, 1)
        t = P(Fc(i, :), :);
        pts = zeros(0, 2);
        for j = 1:3
            p = t(j, :);  q = t(mod(j, 3) + 1, :);
            if (p(1) - x0) * (q(1) - x0) < 0
                s = (x0 - p(1)) / (q(1) - p(1));
                w = p + s * (q - p);
                pts(end+1, :) = w(2:3);         %#ok<AGROW>
            end
        end
        if size(pts, 1) == 2
            segs(end+1, :, :) = pts;            %#ok<AGROW>
        end
    end
end

% =========================================================================
function poly = ordenar_angular(segs, c)
%ORDENAR_ANGULAR Contorno ordenado por angulo alrededor de c.
% Valido solo para contornos estrellados respecto de c (es el caso de la leva).
    pts = unique(round(reshape(permute(segs, [1 2 3]), [], 2), 6), 'rows');
    [~, k] = sort(atan2(pts(:,2) - c(2), pts(:,1) - c(1)));
    poly = pts(k, :);
end

% =========================================================================
function r = radio_por_rayo(segs, c, alphas)
%RADIO_POR_RAYO Radio del contorno para cada direccion desde c.
% Trabaja sobre los SEGMENTOS crudos del corte, sin encadenarlos: asi no hay
% forma de introducir cuerdas falsas. Requiere que el contorno sea estrellado
% respecto de c, cosa que se verifica al final (r sin huecos).
    r = zeros(numel(alphas), 1);
    A = squeeze(segs(:, 1, :)) - c;      % (n,2) extremo inicial
    B = squeeze(segs(:, 2, :)) - c;      % (n,2) extremo final
    for k = 1:numel(alphas)
        d = [cos(alphas(k)) sin(alphas(k))];
        m = [-d(2) d(1)];
        sa = A * m';  sb = B * m';
        hit = (sa <= 0 & sb > 0) | (sb <= 0 & sa > 0);
        ia = find(hit);
        rr = -inf;
        for j = ia'
            t = sa(j) / (sa(j) - sb(j));
            w = A(j, :) + t * (B(j, :) - A(j, :));
            rr = max(rr, w * d');
        end
        r(k) = rr;
    end
end

% =========================================================================
function informe(G)
    fprintf('\n=== geometria medida del CAD (%s) ===\n', G.assets);
    fprintf('  pivote O        = (%.3f, %.3f) mm   eje = [%g %g %g]\n', G.O, G.eje);
    fprintf('  eje de leva C   = (%.3f, %.3f) mm\n', G.C);
    fprintf('\n  --- masas (densidad 1, escalar por rho) ---\n');
    for k = {'balancin','leva','sop_der','sop_izq'}
        s = G.(k{1});
        fprintf('  %-9s V = %9.1f mm^3   centroide = (%7.3f %7.3f %7.3f)\n', ...
                k{1}, s.V_mm3, s.C_mm);
    end
    fprintf('\n  --- leva ---\n');
    fprintf('  disco en x = [%.2f %.2f] mm (espesor %.2f)\n', ...
            G.leva.x_disco, diff(G.leva.x_disco));
    fprintf('  perfil: %d puntos, r = %.3f .. %.3f mm respecto de C\n', ...
            size(G.leva.poligono,1), G.leva.r_min, G.leva.r_max);
    fprintf('\n  --- balancin ---\n');
    fprintf('  cara del seguidor en y = %.3f mm; punta del brazo en z = %.3f mm\n', ...
            G.balancin.d_cara, G.balancin.z_punta);
    fprintf('  mazo: y = %.2f mm, base en z = %.2f mm\n', ...
            G.balancin.y_mazo, G.balancin.z_mazo);
    fprintf('  region de contacto: %d puntos\n', size(G.balancin.contacto,1));
    fprintf('\n  --- cinematica ---\n');
    fprintf('  rho(theta) = %.0f cos(th) + %.0f sin(th) - %.1f  [mm]\n', ...
            G.C(1), G.C(2), G.balancin.d_cara);
    fprintf('  rho maximo alcanzable = %.3f mm (en theta = %.2f deg)\n', ...
            G.rho_max_alcanzable, rad2deg(G.theta_rho_max));
    if G.leva.r_max > G.rho_max_alcanzable
        fprintf('  AVISO: la nariz (r=%.2f) excede en %.2f mm lo que la cara alcanza\n', ...
                G.leva.r_max, G.leva.r_max - G.rho_max_alcanzable);
        fprintf('         -> la suelta ocurre por resbalon en la punta del brazo\n');
    end
    fprintf('  theta con la leva en su circulo base: %.2f deg (holgura de suelta)\n', ...
            rad2deg(theta_de_rho(G, G.leva.r_min)));
end

function th = theta_de_rho(G, r)
    amp = hypot(G.C(1), G.C(2));
    v = (r + G.balancin.d_cara) / amp;
    v = max(min(v, 1), -1);
    th = G.theta_rho_max - acos(v);
end
