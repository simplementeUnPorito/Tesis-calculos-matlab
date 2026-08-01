# Física loca — martillo de leva en Simscape Multibody

Modelo dinámico del martillo de leva diseñado en `Assets/Fusion 360`, para
dimensionar el mecanismo que genera la fuente sísmica del ensayo MASW.

## Uso

```matlab
cd('C:\Github\Tesis\src\calculos_modelados\matlab\Fisica loca')

extraer_geometria_cad();          % mide el CAD -> geometria_cad.mat
p = params_martillo();            % parámetros físicos
construir_martillo_leva(p);       % genera martillo_leva.slx
R = simular_martillo();           % corre, grafica y reporta
```

Variantes útiles:

```matlab
R = simular_martillo('rpm', 60);                  % otra velocidad de leva
R = simular_martillo('material', 'aluminio');     % otro material
R = simular_martillo('escala', 10);               % lleva el CAD 1:10 a campo
R = simular_martillo('recorte_nariz', 41);        % evalúa recortar la nariz
```

## Archivos

| archivo | qué hace |
|---|---|
| `extraer_geometria_cad.m` | mide ejes, perfil de leva y propiedades másicas desde los STL |
| `params_martillo.m` | material, rpm, rigideces de contacto, escala |
| `construir_martillo_leva.m` | arma `martillo_leva.slx` desde cero (bloques + cableado) |
| `contacto_leva_seguidor.m` | contacto unilateral analítico leva–seguidor |
| `simular_martillo.m` | corre, calcula reacciones y par, grafica y reporta |
| `geometria_cad.mat` | geometría medida (regenerable) |
| `contacto_leva_datos.mat` | sólo arreglos numéricos, para `coder.load` en el bloque |

Los `.mat` y el `.slx` son regenerables; se pueden borrar sin perder nada.

## Geometría medida del CAD

El CAD está pensado para imprimir en 3D, con escala de maqueta pero formas
definitivas. Todo lo siguiente se **midió** de los archivos, no se asumió:

- El mecanismo es plano y vive en el plano **YZ**; **z es la vertical** y los
  dos ejes de giro son **+X** (los pernos).
- **Pivote O = (y,z) = (0, 0)**. Verificado por dos caminos independientes: el
  buje ø10 de `Hammer.step` está en el origen, y `soporte1.step` tiene su
  agujero de perno ø5 en `[46, 0, 0]`.
- **Eje de la leva C = (y,z) = (25, 40)**. Igual: buje de `Leva.step` más el
  agujero de `soporte1.step` en `[46, 25, 40]`.
- `C = (25,40)` es exactamente **1:10** de `Lt=250, Ht=400` del modelo Python
  `src/calculos_modelados/python/leva_martillo/leva_martillo.py`.
- **θ = 0 es el instante del impacto**: en esa pose la base del mazo queda en
  `z = −15 mm`, el mismo nivel que la base de los soportes (el piso). θ > 0
  levanta el mazo. El ensamblaje STEP está dibujado en **θ = +15.897°**, que
  es una pose ya levantada (sale de la transformada de la ocurrencia
  `Hammer:1`, una rotación pura alrededor de +X sin componente y/z, lo que
  confirma que el eje de pivote pasa por y=z=0).
- Balancín en L: brazo del mazo ~95 mm sobre +y con la cabeza en T; brazo
  seguidor vertical de 61.6 mm, con la **cara plana en y = +5 mm**.
- Leva tipo caracol: disco de **10 mm** de espesor (x ∈ [−5, +5]; el resto del
  ancho es el eje), radio **17.125 a 45.108 mm** respecto de C.
- Distancia del eje de la leva a la cara del seguidor:
  **ρ(θ) = 25·cos θ + 40·sin θ − 5** [mm], la misma forma que usa el modelo
  Python, con el término −5 por el desplazamiento de la cara.

Masas con PLA (1240 kg/m³): balancín 21.85 g, leva 18.47 g, soportes 44.45 g
cada uno.

## Decisiones de modelado y por qué

**Por qué no `smimport`.** En R2025b `smimport` sólo acepta XML de Simscape
Multibody Link o URDF; no lee STEP. Y Fusion 360 no tiene plugin de Multibody
Link (existe para SolidWorks, Inventor y Creo). Por eso el modelo se construye
bloque por bloque con `construir_martillo_leva.m`.

**Por qué STL y no STEP en los `File Solid`.** El bloque no logra abrir los
`.step` en esta instalación (reporta *file does not exist* con cualquier forma
de la ruta) y sí abre los `.stl`. Como además las propiedades másicas se miden
de esas mismas mallas, el modelo queda auto-consistente.

**Por qué el contacto es analítico y no `Spatial Contact Force`.** El perfil de
la leva **no es convexo**: se desvía 3.66 mm de su envolvente convexa
justamente en el escalón que produce la suelta. `Spatial Contact Force`
convexifica la geometría, así que borraría el escalón y la leva nunca soltaría.
En cambio la leva **sí es estrellada** respecto de C, o sea que su contorno se
describe con un radio univaluado `r(α)`; `contacto_leva_seguidor.m` transforma
cada punto de la región de contacto del seguidor al marco de la leva y compara
radios, con costo O(1) por punto.

**Por qué las juntas casi no usan puertos de sensado.** El orden de los puertos
de una junta de Simscape depende de qué opciones estén habilitadas y no hay
forma de consultar su nombre en tiempo de edición. Para no depender de eso,
sólo se sacan las cuatro señales que necesita el contacto (θ, ω_bal, ψ, ω_cam)
y todo lo demás se obtiene por otras vías.

**Por qué las reacciones y el par se calculan en post-proceso.** El mecanismo es
plano y de 1 grado de libertad, así que las reacciones en O y en C y el par del
motor salen exactos de la dinámica del cuerpo rígido a partir de θ, ω, α y la
fuerza de contacto. Es equivalente a sensarlas y no depende del orden de
puertos.

**Accionamiento de la leva.** Prescribir la posición de una junta exige dos
derivadas del input y el bloque `Simulink-PS Converter` sólo admite una
explícita, así que se usa su filtro de 2º orden (`p.tau_filtro = 0.5 ms`). Para
una rampa de velocidad constante es exacto en régimen; sólo deja un retardo
fijo de `w_cam · tau` en ψ (≈0.12° a 40 rpm).

**El piso está modelado como límite inferior de la junta**, no como contacto
explícito con el mazo. Consecuencia importante: mientras el balancín está
apoyado contra el tope, la reacción calculada en O absorbe también la fuerza del
piso y **no representa la carga del cojinete**. `simular_martillo` excluye esas
muestras del máximo que reporta y lo avisa. Para dimensionar el eje en el
instante del golpe hay que modelar el yunque como contacto explícito.

## Resultado principal: la leva y el balancín son coplanares y chocan

**El plato del balancín y el disco de la leva ocupan exactamente el mismo rango
axial `x ∈ [−5, +5] mm`.** Medido de las mallas, por región:

| pieza / región | extensión en x |
|---|---|
| brazo del mazo (y>30, z<20) | −5.00 … +5.00 |
| cabeza del mazo (y>80) | −5.00 … +5.00 |
| brazo seguidor (z>25) | −5.00 … +5.00 |
| **disco de la leva** (r>20 de C) | **−5.00 … +5.00** |
| ejes / pernos | −56.6 … +56.6 |

Como no hay desfase axial, el cuerpo de la leva barre el espacio que ocupa el
brazo del mazo. El círculo base (r = 17.13 mm) ya alcanza al brazo a partir de
**θ ≈ 28.6°**, y la nariz lo alcanza antes:

| θ [deg] | 0 | 15 | 25 | **30** | 35 | 40 | 45 |
|---|---|---|---|---|---|---|---|
| dist. de C al eje del brazo [mm] | 40.0 | 32.2 | 25.7 | 22.1 | 18.4 | 14.6 | 10.6 |
| ¿solapa con el brazo? | no | no | no | **sí** | sí | sí | sí |

Ver `interferencia_leva_brazo.png`: a θ=45° la leva atraviesa el brazo.

Esto coincide con lo que exige el modelo Python, que pone `z_leva = 0` y
`z_brazo = 30` y avisa que *"proyectados sobre el mismo plano se cruzan
(dist 2D = 48.8 mm): el desfase axial es obligatorio"*.

**Corrección sugerida:** separar axialmente la leva del brazo del mazo, al menos
`(e_leva + e_brazo)/2 = 10 mm` de eje a eje según el criterio del script Python
(que pide ≥16 mm y usa 30 mm). O sea: correr el brazo del mazo a otro plano y
dejar el brazo seguidor compartiendo plano con la leva, que es donde debe haber
contacto.

Mientras la interferencia exista, los números de impacto, par y reacciones que
sale del modelo **no son utilizables**: el contacto detecta la penetración
leva-brazo (hasta 15 mm) y la convierte en golpes espurios. La geometría hay que
arreglarla antes de dimensionar nada.

### Segundo punto a revisar: rango útil del seguidor de cara plana

Con `ρ(θ) = 25·cos θ + 40·sin θ − 5`, el punto de contacto sobre el brazo
seguidor está a `40·cos θ − 25·sin θ` mm del pivote: 40 mm a θ=0, 14 mm a
θ=40°, y **0 a θ=58°**. O sea que al subir, el contacto se desliza hacia el
pivote y el brazo de palanca se anula. ρ tiene su máximo de 42.17 mm justo en
θ=58°, mientras la leva del CAD llega a **45.11 mm** de radio. Conviene revisar
si la nariz no está sobredimensionada respecto del rango útil (el parámetro
`recorte_nariz` está para evaluarlo), aunque esto **no se pudo verificar en
simulación todavía** porque la interferencia anterior domina el resultado.

## Trampas encontradas (para no repetirlas)

- **Encadenar segmentos en un polígono es frágil.** La primera versión de la
  extracción cerraba el contorno de la leva antes de tiempo y lo tapaba con una
  cuerda recta de 42.7 mm (contra una mediana de lado de 0.894 mm). Esa cuerda
  falsa actuaba como un flanco vertical inexistente y en la simulación golpeaba
  al seguidor con penetraciones de 9 mm y ω de 575 rad/s, que se parecía mucho a
  un defecto de diseño real. Ahora `r(α)` se calcula intersecando rayos contra
  los **segmentos crudos**, sin encadenarlos, y `extraer_geometria_cad` valida
  el perfil (huecos y pendiente máxima) antes de dejar simular.
- `coder.load` carga **todas** las variables del `.mat`, así que los datos del
  contacto van en un archivo aparte, sin el struct `G` ni function handles.
- `simlog` guarda por defecto sólo las **últimas 5000 muestras**
  (`SimscapeLogLimitData`), lo que hace parecer que la simulación empieza tarde.
- La licencia de Simscape Multibody de esta instalación responde al nombre de
  feature **legado** `SimMechanics`; `license('test','Simscape_Multibody')`
  devuelve 0 aunque el producto funcione.
- **Elegir un solo punto de contacto (el más penetrado) cuelga el solver.** La
  fuerza salta cada vez que cambia cuál es el punto más profundo entre los
  muestreados, y `daessc` se clava reduciendo el paso indefinidamente (una
  corrida de 3.3 s no terminó en 20 min). La versión actual reparte la rigidez
  entre los N puntos y suma las contribuciones, así cada punto entra y sale del
  contacto de forma continua; con eso la corrida tarda ~1 min.

## Pendiente

- **`gap_brazo` de `contacto_leva_seguidor.m` no es confiable**: devuelve un
  valor constante (−5 mm) en vez de la separación real. El chequeo de
  interferencia de este README se hizo con un script geométrico aparte, no con
  esa salida. Hay que reescribirlo como distancia polígono-polígono.
- Excluir el brazo del mazo de la región de contacto del seguidor una vez que la
  geometría tenga el desfase axial (hoy `contacto` toma todos los puntos con
  `z > 2 & y > −1`, que incluye el techo del brazo del mazo hasta y=95).
- Modelar el yunque como contacto explícito en vez de límite de junta, para
  poder dimensionar los ejes con la carga del golpe.
- Comparar contra `leva_martillo.py` recién cuando la geometría no interfiera.

## Relación con el modelo Python

`leva_martillo.py` resuelve el mismo mecanismo en 2D con el perfil de leva
calculado como envolvente. Dos diferencias que conviene tener presentes:

1. El perfil de la leva del CAD **no** es el que calcula el script: el CAD tiene
   círculo base r = 17.1 mm, mientras el script da 88.1 mm (8.8 mm a escala
   1:10). El CAD se dibujó con otro juego de parámetros.
2. El script reporta **59.25 J de energía de impacto contra 34.99 J de energía
   potencial liberada**, lo que viola conservación de energía; el valor
   coherente con esos 34.99 J sería ω = 3.89 rad/s (v = 3.30 m/s), no los
   4.30 m/s que informa. Conviene revisar `v_imp` en `simular()` antes de usar
   esos números.
