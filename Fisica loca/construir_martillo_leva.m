function mdl = construir_martillo_leva(p, varargin)
%CONSTRUIR_MARTILLO_LEVA Arma el modelo Simscape Multibody del martillo de leva.
%
% Construye martillo_leva.slx desde cero: cuerpos con la geometria real de los
% STEP de Fusion 360, juntas revolutas en los ejes medidos del CAD, la leva
% accionada a velocidad constante, contacto unilateral leva-seguidor y el piso
% como tope inferior del balancin.
%
% El mecanismo es plano en YZ y los dos ejes de giro son +X. Como una junta
% revoluta gira alrededor de SU eje Z, cada junta se orienta con un Rigid
% Transform que lleva Z local a X global (giro de +90 deg alrededor de +Y).
%
% Solo se sacan por puertos de sensado las cuatro senales que necesita el
% contacto (theta, w_bal, psi, w_cam). Todo lo demas -par del motor, reacciones
% en los cojinetes, fuerza del tope- se lee despues de simlog, que identifica
% las variables por nombre y no depende del orden de los puertos.
%
%   mdl = construir_martillo_leva(params_martillo());
%
% Ver tambien: params_martillo, extraer_geometria_cad, simular_martillo

    q = inputParser();
    q.addParameter('Nombre', 'martillo_leva', @(s) ischar(s) || isstring(s));
    q.addParameter('Abrir', false, @islogical);
    q.parse(varargin{:});
    o = q.Results;

    mdl  = char(o.Nombre);
    aqui = fileparts(mfilename('fullpath'));
    G    = p.G;

    if bdIsLoaded(mdl), close_system(mdl, 0); end
    new_system(mdl);
    load_system(mdl);
    B = @(n) [mdl '/' n];

    % =====================================================================
    % infraestructura
    % =====================================================================
    add(B('Solver'),  'nesl_utility/Solver Configuration',        [  60  540]);
    add(B('World'),   'sm_lib/Frames and Transforms/World Frame', [  60  300]);
    add(B('MechCfg'), 'sm_lib/Utilities/Mechanism Configuration', [  60  420]);
    set_param(B('MechCfg'), 'GravityVector', mat2str(p.gravedad), ...
                            'GravityVectorUnits', 'm/s^2');

    % =====================================================================
    % soportes, rigidos al mundo
    % =====================================================================
    add(B('Sop der'), 'sm_lib/Body Elements/File Solid', [ 260  100]);
    add(B('Sop izq'), 'sm_lib/Body Elements/File Solid', [ 260  180]);
    solido(B('Sop der'), fullfile(G.assets, 'Base derecha.stl'),   p.rho_soporte, p, [.75 .75 .78]);
    solido(B('Sop izq'), fullfile(G.assets, 'Base izquierda.stl'), p.rho_soporte, p, [.75 .75 .78]);

    % =====================================================================
    % balancin: O esta en el origen del CAD, no hace falta trasladar
    % =====================================================================
    add(B('EjeO'),     'sm_lib/Frames and Transforms/Rigid Transform', [ 220  300]);
    add(B('Pivote'),   'sm_lib/Joints/Revolute Joint',                 [ 340  300]);
    add(B('EjeO inv'), 'sm_lib/Frames and Transforms/Rigid Transform', [ 480  300]);
    add(B('Balancin'), 'sm_lib/Body Elements/File Solid',              [ 600  300]);
    girar_z_a_x(B('EjeO'),     +90);
    girar_z_a_x(B('EjeO inv'), -90);
    solido(B('Balancin'), fullfile(G.assets, 'Mazo.stl'), p.rho_balancin, p, [.30 .45 .75]);

    % theta = 0 es el impacto: ahi la base del mazo llega al nivel del piso.
    % El piso/yunque se modela como el limite inferior de la junta.
    set_param(B('Pivote'), ...
        'LowerLimitSpecify', 'on', ...
        'LowerLimitBound', '0', 'LowerLimitBoundUnits', 'deg', ...
        'LowerLimitStiffness', 'p.tope_rigidez', ...
        'LowerLimitStiffnessUnits', 'N*m/deg', ...
        'LowerLimitDamping', 'p.tope_amort', ...
        'LowerLimitDampingUnits', 'N*m/(deg/s)', ...
        'LowerLimitTransitionRegionWidth', 'p.tope_ancho', ...
        'LowerLimitTransitionRegionWidthUnits', 'deg', ...
        'UpperLimitSpecify', 'on', ...
        'UpperLimitBound', 'p.tope_sup_deg', 'UpperLimitBoundUnits', 'deg', ...
        'UpperLimitStiffness', 'p.tope_rigidez', ...
        'UpperLimitStiffnessUnits', 'N*m/deg', ...
        'UpperLimitDamping', 'p.tope_amort', ...
        'UpperLimitDampingUnits', 'N*m/(deg/s)', ...
        'UpperLimitTransitionRegionWidth', 'p.tope_ancho', ...
        'UpperLimitTransitionRegionWidthUnits', 'deg', ...
        'DampingCoefficient', 'p.amort_pivote', ...
        'DampingCoefficientUnits', 'N*m/(rad/s)', ...
        'TorqueActuationMode', 'InputTorque', ...
        'SensePosition', 'on', 'SenseVelocity', 'on', ...
        'PositionTargetSpecify', 'on', ...
        'PositionTargetValue', 'p.theta_ini_deg', ...
        'PositionTargetValueUnits', 'deg', ...
        'PositionTargetPriority', 'High');

    % =====================================================================
    % leva: traslada a C, gira Z->X, junta, y deshace para colocar el solido
    % =====================================================================
    add(B('EjeC'),     'sm_lib/Frames and Transforms/Rigid Transform', [ 220  700]);
    add(B('Leva j'),   'sm_lib/Joints/Revolute Joint',                 [ 340  700]);
    add(B('EjeC inv'), 'sm_lib/Frames and Transforms/Rigid Transform', [ 480  700]);
    add(B('Leva'),     'sm_lib/Body Elements/File Solid',              [ 600  700]);

    C_m = [0 G.C(1) G.C(2)] * 1e-3 * p.escala;
    girar_z_a_x(B('EjeC'), +90);
    set_param(B('EjeC'), 'TranslationMethod', 'Cartesian', ...
        'TranslationCartesianOffset', mat2str(C_m), ...
        'TranslationCartesianOffsetUnits', 'm');
    % el inverso: primero deshace el giro, luego destraslada en coords del CAD
    girar_z_a_x(B('EjeC inv'), -90);
    set_param(B('EjeC inv'), 'TranslationMethod', 'Cartesian', ...
        'TranslationCartesianOffset', mat2str(-C_m), ...
        'TranslationCartesianOffsetUnits', 'm');
    solido(B('Leva'), fullfile(G.assets, 'Leva.stl'), p.rho_leva, p, [.85 .55 .20]);

    set_param(B('Leva j'), ...
        'MotionActuationMode', 'InputMotion', ...
        'TorqueActuationMode', 'ComputedTorque', ...
        'DampingCoefficient', 'p.amort_leva', ...
        'DampingCoefficientUnits', 'N*m/(rad/s)', ...
        'SensePosition', 'on', 'SenseVelocity', 'on');

    % =====================================================================
    % accionamiento: psi = w*t, con la derivada exacta
    % =====================================================================
    add(B('psi rampa'), 'simulink/Sources/Ramp',              [  60  760]);
    add(B('psi->PS'),   'nesl_utility/Simulink-PS Converter', [ 200  790]);
    set_param(B('psi rampa'), 'slope', 'p.w_cam', 'start', '0', 'InitialOutput', '0');
    % Prescribir la posicion de una junta exige dos derivadas del input, y este
    % bloque solo admite UNA explicita. Se usa entonces el filtro de 2do orden,
    % que las genera internamente: para una rampa de velocidad constante es
    % exacto en regimen y solo deja un retardo fijo de w_cam*tau en psi
    % (con tau = 0.5 ms y 40 rpm son ~0.12 deg).
    set_param(B('psi->PS'), 'FilteringAndDerivatives', 'filter', ...
        'SimscapeFilterOrder', '2', ...
        'InputFilterTimeConstant', 'p.tau_filtro', 'Unit', 'rad');

    % =====================================================================
    % contacto leva-seguidor
    % =====================================================================
    add(B('th->SL'), 'nesl_utility/PS-Simulink Converter', [ 700  250]);
    add(B('wb->SL'), 'nesl_utility/PS-Simulink Converter', [ 700  310]);
    add(B('ps->SL'), 'nesl_utility/PS-Simulink Converter', [ 700  650]);
    add(B('wc->SL'), 'nesl_utility/PS-Simulink Converter', [ 700  710]);
    set_param(B('th->SL'), 'Unit', 'rad');    set_param(B('wb->SL'), 'Unit', 'rad/s');
    set_param(B('ps->SL'), 'Unit', 'rad');    set_param(B('wc->SL'), 'Unit', 'rad/s');

    add(B('k'),      'simulink/Sources/Constant', [ 700  390]);
    add(B('c'),      'simulink/Sources/Constant', [ 700  440]);
    add(B('escala'), 'simulink/Sources/Constant', [ 700  490]);
    add(B('r max'),  'simulink/Sources/Constant', [ 700  540]);
    set_param(B('k'),      'Value', 'p.k_contacto');
    set_param(B('c'),      'Value', 'p.c_contacto');
    set_param(B('escala'), 'Value', 'p.escala');
    set_param(B('r max'),  'Value', 'p.r_leva_max');

    add(B('Contacto'), 'simulink/User-Defined Functions/MATLAB Function', [ 880  300]);
    set_param(B('Contacto'), 'Position', [880 300 1040 480]);
    script_contacto(B('Contacto'));

    add(B('F bal'),  'sm_lib/Forces and Torques/External Force and Torque', [1300  260]);
    add(B('F leva'), 'sm_lib/Forces and Torques/External Force and Torque', [1300  660]);
    for b = {B('F bal'), B('F leva')}
        set_param(b{1}, 'EnableForceY', 'on', 'EnableForceZ', 'on', ...
                        'EnableTorqueX', 'on', ...
                        'ForceResolutionFrame', 'World', ...
                        'TorqueResolutionFrame', 'World');
    end

    % Simulink -> physical signal. El orden de LConn de External Force and
    % Torque sigue el orden del dialogo: fy, fz, tx.
    ps = {'Fy','N',   'Fz','N',   'Tb','N*m', ...
          'Fy2','N',  'Fz2','N',  'Tc','N*m'};
    for i = 1:2:numel(ps)
        n = ['->PS ' ps{i}];
        add(B(n), 'nesl_utility/Simulink-PS Converter', [1150  200 + 60*((i+1)/2)]);
        set_param(B(n), 'Unit', ps{i+1});
    end

    % =====================================================================
    % registro de las senales del contacto
    % =====================================================================
    add(B('Diag'), 'simulink/Signal Routing/Mux', [1150  620]);
    set_param(B('Diag'), 'Inputs', '4', 'Position', [1150 600 1155 680]);
    add(B('log'), 'simulink/Sinks/To Workspace', [1250  620]);
    set_param(B('log'), 'VariableName', 'diag', ...
        'SaveFormat', 'Structure With Time', 'SampleTime', '-1');

    % =====================================================================
    % cableado
    % =====================================================================
    L = @(b, i) get_param(B(b), 'PortHandles').LConn(i);
    R = @(b, i) get_param(B(b), 'PortHandles').RConn(i);
    In = @(b, i) get_param(B(b), 'PortHandles').Inport(i);
    Out = @(b, i) get_param(B(b), 'PortHandles').Outport(i);
    une = @(a, b) add_line(mdl, a, b, 'autorouting', 'smart');

    % red fisica: el mundo alimenta las dos cadenas, los soportes,
    % la configuracion del mecanismo y el solver
    une(R('World',1), L('EjeO',1));
    une(R('World',1), L('EjeC',1));
    une(R('World',1), R('Sop der',1));
    une(R('World',1), R('Sop izq',1));
    une(R('World',1), R('MechCfg',1));
    une(R('World',1), R('Solver',1));

    % cadena del balancin:  B=LConn1, t=LConn2, F=RConn1, q=RConn2, w=RConn3
    une(R('EjeO',1),     L('Pivote',1));
    une(R('Pivote',1),   L('EjeO inv',1));
    une(R('EjeO inv',1), R('Balancin',1));

    % cadena de la leva:  B=LConn1, movimiento=LConn2
    une(R('EjeC',1),     L('Leva j',1));
    une(R('Leva j',1),   L('EjeC inv',1));
    une(R('EjeC inv',1), R('Leva',1));

    % accionamiento
    une(Out('psi rampa',1), In('psi->PS',1));
    une(R('psi->PS',1),     L('Leva j',2));

    % sensado hacia el contacto
    une(R('Pivote',2), L('th->SL',1));
    une(R('Pivote',3), L('wb->SL',1));
    une(R('Leva j',2), L('ps->SL',1));
    une(R('Leva j',3), L('wc->SL',1));

    une(Out('th->SL',1), In('Contacto',1));
    une(Out('wb->SL',1), In('Contacto',2));
    une(Out('ps->SL',1), In('Contacto',3));
    une(Out('wc->SL',1), In('Contacto',4));
    une(Out('k',1),      In('Contacto',5));
    une(Out('c',1),      In('Contacto',6));
    une(Out('escala',1), In('Contacto',7));
    une(Out('r max',1),  In('Contacto',8));

    % fuerzas de contacto hacia los cuerpos
    nombres = {'Fy','Fz','Tb','Fy2','Fz2','Tc'};
    for i = 1:6
        une(Out('Contacto', i), In(['->PS ' nombres{i}], 1));
    end
    une(R('->PS Fy',1),  L('F bal',1));
    une(R('->PS Fz',1),  L('F bal',2));
    une(R('->PS Tb',1),  L('F bal',3));
    une(R('->PS Fy2',1), L('F leva',1));
    une(R('->PS Fz2',1), L('F leva',2));
    une(R('->PS Tc',1),  L('F leva',3));

    % la fuerza sobre el balancin se aplica en O (marco F del pivote) y la
    % reaccion sobre la leva en C (marco F de la junta de la leva)
    une(R('F bal',1),  R('Pivote',1));
    une(R('F leva',1), R('Leva j',1));

    % diagnostico: delta, gap, py, pz
    for i = 1:4
        une(Out('Contacto', 6 + i), In('Diag', i));
    end
    une(Out('Diag',1), In('log',1));

    % =====================================================================
    % solver
    % =====================================================================
    set_param(mdl, ...
        'SolverType', 'Variable-step', 'Solver', 'daessc', ...
        'RelTol', num2str(p.tol_rel), 'MaxStep', num2str(p.dt_max), ...
        'StopTime', num2str(p.t_fin), ...
        'SimscapeLogType', 'all', 'SimscapeLogName', 'simlog', ...
        'ReturnWorkspaceOutputs', 'on');

    destino = fullfile(aqui, [mdl '.slx']);
    save_system(mdl, destino);
    fprintf('-> %s\n', destino);
    if o.Abrir, open_system(mdl); end
end

% =========================================================================
function add(dest, src, pos)
    add_block(src, dest, 'MakeNameUnique', 'off');
    pp = get_param(dest, 'Position');
    w = pp(3) - pp(1);  h = pp(4) - pp(2);
    set_param(dest, 'Position', [pos(1) pos(2) pos(1)+w pos(2)+h]);
end

function solido(blk, archivo, rho, p, color)
%SOLIDO Configura un File Solid: geometria de archivo, inercia por densidad.
    if ~isfile(archivo)
        error('construir_martillo_leva:SinCAD', 'Falta el archivo CAD: %s', archivo);
    end
    % El CAD esta en mm. La escala CAD -> maquina real se aplica leyendo el
    % archivo en otra unidad, que es lo unico que ofrece el bloque File Solid.
    switch p.escala
        case 1,    u = 'mm';
        case 10,   u = 'cm';
        case 100,  u = 'dm';
        case 1000, u = 'm';
        otherwise
            error('construir_martillo_leva:Escala', ...
                ['escala = %g no se puede aplicar a la geometria: File Solid solo ' ...
                 'admite unidades de archivo, asi que la escala debe ser 1, 10, ' ...
                 '100 o 1000 (mm, cm, dm, m).'], p.escala);
    end
    set_param(blk, ...
        'ExtGeomFileName', archivo, ...
        'UnitType', 'Custom', 'ExtGeomFileUnits', u, ...
        'InertiaType', 'CalculateFromGeometry', ...
        'BasedOnType', 'Density', ...
        'Density', num2str(rho), 'DensityUnits', 'kg/m^3', ...
        'GraphicVisPropType', 'SimpleVisualProperties', ...
        'GraphicDiffuseColor', mat2str(color));
end

function girar_z_a_x(blk, ang)
%GIRAR_Z_A_X Rigid Transform que lleva el eje Z local al X global (+90 en +Y).
    set_param(blk, 'RotationMethod', 'StandardAxis', ...
                   'RotationStandardAxis', '+Y', ...
                   'RotationAngle', num2str(ang), ...
                   'RotationAngleUnits', 'deg');
end

function script_contacto(blk)
%SCRIPT_CONTACTO Cuerpo del bloque MATLAB Function que envuelve el contacto.
    txt = [ ...
      "function [Fy, Fz, Tb, Fy2, Fz2, Tc, delta, gap, py, pz] = ..."
      "        fcn(theta, w_bal, psi, w_cam, k, c, escala, r_max)"
      "%#codegen"
      "% La reaccion sobre la leva es igual y opuesta a la fuerza sobre el balancin."
      "[Fy, Fz, Tb, Tc, delta, gap, py, pz] = ..."
      "    contacto_leva_seguidor(theta, psi, w_bal, w_cam, k, c, escala, r_max);"
      "Fy2 = -Fy;"
      "Fz2 = -Fz;"
    ];
    S = sfroot();
    ch = S.find('-isa', 'Stateflow.EMChart', 'Path', blk);
    ch.Script = char(strjoin(txt, newline));
end
