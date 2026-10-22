"""Modelo termico de la placa BUCK LM2596 sobre la geometria STEP importada.

Parte de BUCK_LM2596_termico.mph (geometria importada y arreglada en COMSOL) y
guarda TODO en un archivo nuevo, BUCK_LM2596_termico_STEP.mph. El original no se
modifica.

Pasos:
  1. Geometria: se conservan la placa y los componentes del STEP. Las pistas, pads,
     zonas de cobre y barriles de via del STEP (laminas de 25-40 um) se quitan, porque
     mallarlas no cabe en memoria (la malla fina llego a 11.5 GB y el sistema la mato).
     En su lugar:
       - el cobre de F.Cu y B.Cu entra como capa delgada de 35 um con la forma real
         (mascaras exportadas de KiCad, geometria/mascara_cobre_*.txt);
       - cada via es un cilindro solido con la conductividad equivalente de su barril;
       - soldadura (SAC305) entre la placa y el tab de U1, las terminales de D1 y los
         pines de U1 (los modelos 3D de KiCad empiezan en z = 1.595 mm);
       - pads de L1 y L2 como bloques de cobre (sus cuerpos no se modelan).
  2. Materiales por forma y tamano de cada dominio.
  3. Heat Transfer in Solids: fuentes como "Heat rate" (U1 en la cara tab-cuerpo,
     D1 en sus terminales, L1 y L2 en sus pads: los cuerpos de L1/L2 no se modelan),
     conveccion natural con correlaciones explicitas segun la orientacion de la cara,
     radiacion opcional.
  4. Malla: tetraedros con la geometria escalada en z (soldadura y pads delgados).
  5. Estudios: caso nominal, barrido de corriente y variantes.
Coordenadas: las del STEP de KiCad (mm), x 20..80, y -50..-20.
"""
from __future__ import annotations

import csv
import math
import os
from pathlib import Path

import mph
import numpy as np

FEM = Path(__file__).resolve().parent.parent
ORIGEN = FEM / 'BUCK_LM2596_termico.mph'
# MODO=capa (por defecto): pistas, pads y plano como capa delgada (malla ligera).
# MODO=cobre: pistas, pads y plano GND como solidos del STEP (mas pesado de mallar).
# MODO=aire: como "capa", pero con una caja de aire alrededor y el flujo de convección
#            natural calculado (laminar + gravedad + Non-Isothermal Flow), en el caso
#            12 V -> 2 V, 0.5 A.
MODO = os.environ.get('MODO', 'capa')
DESTINO = FEM / {'cobre': 'BUCK_LM2596_termico_STEP_cobre.mph',
                 'aire': 'BUCK_LM2596_termico_12V_2V_05A_aire.mph'}.get(MODO, 'BUCK_LM2596_termico_STEP.mph')
# Caja de aire (mm): placa en x 20..80, y -50..-20, z 0..1.6; margen extra arriba (pluma)
AIRE = ((0, 100), (-70, 0), (-25, 50))  # 100 x 70 x 75 mm (la caja de 120 x 90 x 110 convergia muy lento)
# Agujeros THT (coordenadas del STEP, diametro en mm). Con aire se rellenan de soldadura,
# como en la placa real; si no, queda un anillo de aire finisimo alrededor de cada pata.
AGUJEROS_THT = [(25.4, -30.55, 0.8), (25.4, -34.05, 0.8), (45.62, -24.6, 0.8), (48.16, -24.6, 0.8),
                (50.7, -24.6, 0.8), (55.5, -29.1, 1.0), (60.5, -29.1, 1.0), (66.4, -29.3, 0.8),
                (69.9, -29.3, 0.8), (22.3, -47.8, 1.5), (22.3, -22.3, 1.5), (77.7, -22.3, 1.5),
                (77.8, -47.8, 1.5)]
FIG = FEM / 'figuras'

GEO = FEM / 'geometria'
Z_PCB, Z_PAD, Z_COMP = 1.51, 1.55, 1.595  # tope del FR4 / tope de un pad / fondo de los componentes
ESCALA_Z = 10  # estiramiento en z para mallar
# Pads de L1 y L2 (coordenadas del STEP)
PADS_L = {
    'pad_l1a': ('Pad L1 1', (42.10, 45.00), (-44.80, -39.40)),
    'pad_l1b': ('Pad L1 2', (52.00, 54.90), (-44.80, -39.40)),
    'pad_l2a': ('Pad L2 1', (58.35, 61.25), (-44.60, -39.20)),
    'pad_l2b': ('Pad L2 2', (68.25, 71.15), (-44.60, -39.20)),
}
# Vias (coordenadas de KiCad); en el STEP y = -y_KiCad
VIAS_KICAD = [(35.95, 28.8), (31.65, 22.4), (35.9, 22.4), (31.6, 25.3), (35.85, 25.3), (40.2, 28.8),
              (39.2, 32.2), (48.7, 30.0), (40.1, 25.3), (35.8, 32.3), (31.7, 28.8), (26.6, 42.8),
              (75.6, 47.8), (24.7, 47.8), (40.15, 22.4), (25.4, 35.7)]
SOLDADURA = {
    'sold_tab': ('Soldadura tab U1', (30.8, 40.8), (-29.23, -21.67)),
    'sold_d1a': ('Soldadura D1 terminal 1', (28.1, 29.95), (-44.45, -41.15)),
    'sold_d1b': ('Soldadura D1 terminal 2', (34.25, 36.1), (-44.45, -41.15)),
}
# Pies de los 5 pines de U1 (x de cada pin, y del pie)
for i, (x0, x1) in enumerate(((31.88, 32.92), (33.58, 34.62), (35.28, 36.32), (36.98, 38.02), (38.68, 39.72)), 1):
    SOLDADURA[f'sold_p{i}'] = (f'Soldadura pin {i} U1', (x0, x1), (-36.6, -35.6))

PARAMETROS = [
    ('Vin', '35[V]', 'Voltaje de entrada (esquematico: 35 V -> 1.2 V)'),
    ('Vout', '1.23[V]', 'Voltaje de salida (minimo del LM2596-ADJ)'),
    ('Iout', '3[A]', 'Corriente de carga (peor caso: maximo del LM2596)'),
    ('fsw', '150[kHz]', 'Frecuencia de conmutacion del LM2596 (datasheet)'),
    ('Vsat', '1.16[V]', 'Caida del switch interno a 3 A (datasheet, tipico)'),
    ('Iq', '5[mA]', 'Corriente de reposo del LM2596 (datasheet, tipico)'),
    ('t_sw', '100[ns]', 'Tiempo de conmutacion subida+bajada (SUPUESTO)'),
    ('Vf', '0.5[V]', 'Caida del diodo SS34 a 3 A (datasheet, maximo)'),
    ('DCR1', '40[mohm]', 'Resistencia del devanado de L1 (SUPUESTO)'),
    ('DCR2', '25[mohm]', 'Resistencia del devanado de L2 (SUPUESTO)'),
    ('incluir_L', '1', '1 = incluir perdidas de L1 y L2, 0 = ignorarlas'),
    ('con_vias', '1', '1 = vias de cobre, 0 = vias como FR4'),
    ('eps_rad', '0', 'Emisividad (0 = solo conveccion, como en la LM317)'),
    ('T_amb', '25[degC]', 'Temperatura ambiente'),
    ('Dcic', '(Vout+Vf)/(Vin-Vsat+Vf)', 'Ciclo de trabajo'),
    ('P_U1', 'Vsat*Iout*Dcic+0.5*Vin*Iout*t_sw*fsw+Vin*Iq', 'Perdidas del LM2596'),
    ('P_D1', 'Vf*Iout*(1-Dcic)', 'Perdidas del diodo SS34'),
    ('P_L1', 'incluir_L*Iout^2*DCR1', 'Perdidas de L1'),
    ('P_L2', 'incluir_L*Iout^2*DCR2', 'Perdidas de L2'),
    ('P_tot', 'P_U1+P_D1+P_L1+P_L2', 'Perdidas totales'),
    ('eficiencia', 'Vout*Iout/(Vout*Iout+P_tot)', 'Eficiencia estimada'),
    ('k_cu', '400[W/(m*K)]', 'Conductividad del cobre'),
    ('k_fr4', '0.3[W/(m*K)]', 'Conductividad del FR4'),
    ('t_cu', '35[um]', 'Espesor del cobre (1 oz)'),
    ('k_via', '122[W/(m*K)]', 'k equivalente de una via: barril de 25 um en un taladro de 0.3 mm'),
    ('k_air', '0.028[W/(m*K)]', 'Conductividad del aire (~50 degC)'),
    ('nu_air', '1.8e-5[m^2/s]', 'Viscosidad cinematica del aire (~50 degC)'),
    ('alpha_air', '2.55e-5[m^2/s]', 'Difusividad termica del aire (~50 degC)'),
    ('Pr_air', '0.71', 'Numero de Prandtl del aire'),
    ('L_h', '10[mm]', 'Longitud caracteristica de caras horizontales (area/perimetro de la placa)'),
    ('L_v', '5[mm]', 'Longitud caracteristica de caras verticales (altura tipica)'),
]

MATERIALES = [
    ('mat_fr4', 'FR4 (placa)', 'k_fr4', '1900', '1369'),
    ('mat_abs', 'Relleno generico (C1-C3, R1)', '0.2', '1050', '1400'),
    ('mat_emc', 'EMC (cuerpo de U1 y D1)', '0.7', '1900', '900'),
    ('mat_cu', 'Cobre (tab y pines de U1, pads de L1/L2)', 'k_cu', '8960', '385'),
    ('mat_via', 'Vias: cobre equivalente (FR4 si con_vias = 0)', 'if(con_vias,k_via,k_fr4)', '8960', '385'),
    ('mat_sold', 'Soldadura SAC305', '58', '7400', '230'),
    # Gas ideal a 1 atm: rho = p M / (R T) = 353.2 kg K/m^3 / T
    ('mat_aire_cfd', 'Aire (gas ideal)', '0.026', '353.2[kg*K/m^3]/T', '1005'),
]


def limpiar(j, comp) -> None:
    """Quita la fisica, estudios y resultados del modelo anterior (selecciones rotas)."""
    for t in [str(x) for x in j.result().tags()]:
        j.result().remove(t)
    for t in [str(x) for x in j.result().dataset().tags()]:
        j.result().dataset().remove(t)
    for t in [str(x) for x in j.study().tags()]:
        j.study().remove(t)
    for t in [str(x) for x in j.sol().tags()]:
        j.sol().remove(t)
    for t in [str(x) for x in comp.physics().tags()]:
        comp.physics().remove(t)
    for t in [str(x) for x in comp.variable().tags()]:
        comp.variable().remove(t)
    for t in [str(x) for x in comp.cpl().tags()]:
        comp.cpl().remove(t)
    for t in [str(x) for x in comp.material().tags()]:
        comp.material().remove(t)
    for t in [str(x) for x in comp.selection().tags()]:
        comp.selection().remove(t)
    for t in [str(x) for x in j.func().tags()]:
        j.func().remove(t)
    for t in [str(x) for x in comp.mesh().tags()]:
        comp.mesh().remove(t)


def preparar_geometria(g) -> None:
    """Quita las laminas de cobre y los barriles de via del STEP y agrega soldadura,
    pads de L1/L2 y vias solidas. Todo se inserta antes de Form Union."""
    tags = [str(t) for t in g.feature().tags()]
    pos = tags.index('fin')
    g.run(tags[pos - 1])
    nombres = [str(n) for n in g.objectNames()]
    claves = ('_via(',) if MODO == 'cobre' else ('_copper(', '_pad(', '_via(')
    quitar = [n for n in nombres if any(k in n for k in claves)]
    for n in nombres if MODO != 'cobre' else []:  # plano de B.Cu: segundo objeto "_PCB" de 35 um
        if '_PCB(' in n:
            bb = list(g.obj(n).getBoundingBox())
            if bb[5] - bb[4] < 0.1:
                quitar.append(n)
    print(f'  Objetos del STEP: {len(nombres)}, se quitan {len(quitar)} (cobre, pads, vias)')
    d = g.create('del_cobre', 'Delete')
    d.label('Quitar laminas de cobre y barriles de via (se modelan como capa delgada)')
    d.selection('input').init()
    d.selection('input').set(quitar)
    g.feature().move('del_cobre', pos)
    pos += 1

    def bloque(tag, etiqueta, x, y, z0, z1):
        nonlocal pos
        b = g.create(tag, 'Block')
        b.label(etiqueta)
        b.set('pos', [f'{x[0]:g}', f'{y[0]:g}', f'{z0:g}'])
        b.set('size', [f'{x[1] - x[0]:g}', f'{y[1] - y[0]:g}', f'{z1 - z0:g}'])
        g.feature().move(tag, pos)
        pos += 1
    z_sold = Z_PAD if MODO == 'cobre' else Z_PCB  # con cobre solido la soldadura va sobre el pad
    for tag, (etiqueta, x, y) in SOLDADURA.items():
        bloque(tag, etiqueta, x, y, z_sold, Z_COMP)
    if MODO != 'cobre':  # con cobre solido los pads de L1/L2 ya vienen en el STEP
        for tag, (etiqueta, x, y) in PADS_L.items():
            bloque(tag, etiqueta, x, y, Z_PCB, Z_PAD)
    if MODO == 'aire':
        # Las patas de los capacitores dentro de los agujeros y bajo la placa generan astillas
        # que el mallador del aire no puede llenar. Se cortan a ras de la cara SUPERIOR: el
        # agujero queda lleno de soldadura y el capacitor se apoya sobre sus munones de pata.
        caps = [n for n in nombres if 'CP_Radial' in n]
        bloque('corte_bajo', 'Bloque para cortar patas bajo la placa', (19, 81), (-51, -19), -2, Z_PCB)
        dif = g.create('cortar_patas', 'Difference')
        dif.label('Cortar patas de capacitores bajo la placa')
        dif.selection('input').set(caps)
        dif.selection('input2').set(['corte_bajo'])
        g.feature().move('cortar_patas', pos)
        pos += 1
        (ax0, ax1), (ay0, ay1), (az0, az1) = AIRE
        bloque('aire', 'Aire alrededor de la placa', (ax0, ax1), (ay0, ay1), az0, az1)
        for i, (xh, yh, dh) in enumerate(AGUJEROS_THT, 1):
            c = g.create(f'tht{i}', 'Cylinder')
            c.label(f'Agujero THT {i} relleno de soldadura')
            c.set('pos', [f'{xh:g}', f'{yh:g}', '0'])
            c.set('r', f'{dh / 2 + 0.25:g}')  # mas ancho que el agujero: lo cubre aunque este corrido
            c.set('h', f'{Z_PCB:g}')
            g.feature().move(f'tht{i}', pos)
            pos += 1
    for i, (xk, yk) in enumerate(VIAS_KICAD, 1):
        c = g.create(f'via{i}', 'Cylinder')
        c.label(f'Via {i}')
        c.set('pos', [f'{xk:g}', f'{-yk:g}', '0'])
        c.set('r', '0.15')
        c.set('h', f'{Z_PCB:g}')
        g.feature().move(f'via{i}', pos)
        pos += 1
    g.run()


def clasificar(g) -> dict[int, tuple[str, str]]:
    mf = g.measureFinal()
    mf.selection().geom(3)
    tabla = {}
    for d in range(1, g.getNDomains() + 1):
        mf.selection().set([d])
        x0, x1, y0, y1, z0, z1 = list(mf.getBoundingBox())
        vol = float(mf.getVolume())
        dx, dy, dz = x1 - x0, y1 - y0, z1 - z0
        cerca = lambda a, b, tol=0.012: abs(a - b) < tol
        if dx > 90:  # la caja de aire mide 100 mm o mas; la placa, 60 mm
            tabla[d] = ('mat_aire_cfd', 'Aire')
        elif dx > 50 and dy > 25 and dz > 1:
            tabla[d] = ('mat_fr4', 'Placa')
        elif (cerca(z0, Z_PCB) or cerca(z0, Z_PAD)) and cerca(z1, Z_COMP):
            tabla[d] = ('mat_sold', 'Soldadura')
        elif cerca(z0, Z_PCB) and cerca(z1, Z_PAD) and MODO != 'cobre':
            tabla[d] = ('mat_cu', 'Pad de L1/L2')
        elif dz <= 0.05:
            tabla[d] = ('mat_cu', 'Lamina de cobre (pista, pad o plano)')
        elif cerca(dz, 4.4, 0.05):
            tabla[d] = ('mat_emc', 'Cuerpo de U1')
        elif cerca(dz, 2.31, 0.05):
            tabla[d] = ('mat_emc', 'Cuerpo de D1')
        elif cerca(z1, 2.865, 0.02):
            tabla[d] = ('mat_cu', 'Tab de U1')
        elif cerca(dz, 2.95, 0.05) or (cerca(dz, 0.5, 0.02) and z0 > 3.5):
            tabla[d] = ('mat_cu', 'Pin de U1')
        elif dz > 9 or cerca(dz, 0.45, 0.02):
            tabla[d] = ('mat_abs', 'Capacitor' if dz > 9 else 'R1')
        elif cerca(dz, Z_PCB, 0.02) and cerca(dx, 0.30, 0.02):
            tabla[d] = ('mat_via', 'Via')
        elif cerca(dz, Z_PCB, 0.02) and MODO == 'cobre':
            tabla[d] = ('mat_cu', f'Barril de agujero THT (diam {dx:.2f} mm)')
        elif cerca(dz, Z_PCB, 0.02) and MODO == 'aire' and 0.75 < dx < 2.3:
            tabla[d] = ('mat_sold', f'Agujero THT relleno de soldadura (diam {dx:.2f} mm)')
        elif cerca(dz, Z_PCB, 0.02) and MODO == 'aire' and dx <= 0.75:
            tabla[d] = ('mat_cu', 'Pata de componente dentro de un agujero')
        elif MODO == 'aire' and z1 < 0.01 and dx < 2:
            tabla[d] = ('mat_cu', 'Pata de componente bajo la placa')
        elif MODO == 'aire' and dx < 15 and dy < 15:
            # Con la caja de aire y los agujeros rellenos, los modelos 3D se parten (cuerpo de
            # capacitor cortado a ras de placa, huecos bajo su tapa): relleno generico
            tabla[d] = ('mat_abs', f'Pieza de componente ({dx:.1f} x {dy:.1f} x {dz:.2f} mm, z0 {z0:.2f})')
            print(f'    dominio {d} como relleno generico: {dx:.2f} x {dy:.2f} x {dz:.3f} mm, z0 = {z0:.3f}')
        else:
            raise ValueError(f'Dominio {d} sin clasificar: {dx:.2f} x {dy:.2f} x {dz:.3f} mm, z0 = {z0:.3f}')
    return tabla


def materiales(comp, tabla) -> None:
    for tag, etiqueta, k, rho, cp in MATERIALES:
        m = comp.material().create(tag, 'Common')
        m.label(etiqueta)
        m.selection().set(sorted(d for d, (t, _) in tabla.items() if t == tag))
        pg = m.propertyGroup('def')
        pg.set('thermalconductivity', [k])
        pg.set('density', [rho])
        pg.set('heatcapacity', [cp])
        if tag == 'mat_aire_cfd':
            pg.set('dynamicviscosity', ['1.85e-5[Pa*s]'])
            pg.set('ratioofspecificheats', ['1.4'])


def caja(comp, tag, etiqueta, x, y, z, dim=2, condicion='intersects'):
    b = comp.selection().create(tag, 'Box')
    b.label(etiqueta)
    b.set('entitydim', str(dim))
    b.set('condition', condicion)
    for k, v in zip(('xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax'), (*x, *y, *z)):
        b.set(k, f'{v:g}')
    return tag


def fisica(comp, tabla, ndom) -> None:
    todos = comp.selection().create('sel_todos', 'Explicit')
    todos.label('Todos los dominios')
    todos.set(list(range(1, ndom + 1)))
    ext = comp.selection().create('sel_ext', 'Adjacent')
    ext.label('Caras exteriores')
    ext.set('input', ['sel_todos'])

    # Fuentes: caras que cruzan una caja delgada en el centro de cada zona
    caja(comp, 'sel_qU1', 'Fuente U1: cara tab-cuerpo', (33, 38), (-28, -24), (2.855, 2.875))
    # D1, L1 y L2: una caja por terminal/pad con la cara de contacto COMPLETA dentro
    # (con 'intersects' tambien entrarian las paredes laterales de soldadura y pads)
    def union(tag, etiqueta, cajas, z):
        partes = []
        for i, (x, y) in enumerate(cajas, 1):
            partes.append(caja(comp, f'{tag}_{i}', f'{etiqueta} {i}', x, y, z, condicion='inside'))
        u = comp.selection().create(tag, 'Union')
        u.label(etiqueta)
        u.set('entitydim', '2')
        u.set('input', partes)
    zc = (Z_COMP - 0.004, Z_COMP + 0.004)
    zp = (Z_PAD - 0.003, Z_PAD + 0.004)
    union('sel_qD1', 'Fuente D1: terminales', [((28.0, 30.05), (-44.55, -41.05)), ((34.15, 36.2), (-44.55, -41.05))], zc)
    union('sel_qL1', 'Fuente L1: pads', [((42.0, 45.1), (-44.9, -39.3)), ((51.9, 55.0), (-44.9, -39.3))], zp)
    union('sel_qL2', 'Fuente L2: pads', [((58.25, 61.35), (-44.7, -39.1)), ((68.15, 71.25), (-44.7, -39.1))], zp)
    for sel in ('sel_qU1', 'sel_qD1', 'sel_qL1', 'sel_qL2'):
        n = len(list(comp.selection(sel).entities()))
        print(f'  {sel}: {n} cara(s)')
        assert n >= 1, f'{sel} quedo vacia'

    ht = comp.physics().create('ht', 'HeatTransfer', 'geom1')
    # Arranque 20 K arriba del ambiente: (T - T_amb)^(1/4) tiene derivada infinita en 0
    ht.feature('init1').set('Tinit', 'T_amb+20[K]')

    # Cobre de F.Cu y B.Cu: capa delgada conductiva de 35 um con la forma real. Las
    # mascaras estan en coordenadas X = x_KiCad - 20, Y = 50 - y_KiCad, y en el STEP
    # y = -y_KiCad, asi que X = x - 20 mm, Y = y + 50 mm.
    caja(comp, 'sel_cara_sup', 'Cara superior de la placa (F.Cu)', (19, 81), (-51, -19), (Z_PCB - 0.004, Z_PCB + 0.004), condicion='inside')
    caja(comp, 'sel_cara_inf', 'Cara inferior de la placa (B.Cu)', (19, 81), (-51, -19), (-0.004, 0.004), condicion='inside')
    capas = (('tl_F', 'Cobre F.Cu (pads y pistas)', 'sel_cara_sup', 'mF'),
             ('tl_B', 'Cobre B.Cu (plano GND)', 'sel_cara_inf', 'mB'))
    for tag, etiqueta, sel, mascara in (capas if MODO != 'cobre' else ()):
        tl = ht.create(tag, 'ThinLayer', 2)
        tl.label(etiqueta)
        tl.selection().named(sel)
        tl.set('LayerType', 'Conductive')
        tl.set('ds', 't_cu')
        tl.set('ks_mat', 'userdef')
        # Ojo: la funcion recibe numeros sin unidad (mm); con [mm] COMSOL los pasa a metros
        k = f'k_fr4+(k_cu-k_fr4)*{mascara}(x/1[mm]-20,y/1[mm]+50)'
        tl.set('ks', [k, '0', '0', '0', k, '0', '0', '0', k])
        for prop, val in (('rhos_mat', 'userdef'), ('rhos', '8960'), ('Cp_s_mat', 'userdef'), ('Cp_s', '385')):
            tl.set(prop, val)
    for tag, etiqueta, sel, pot in (('q_U1', 'U1 LM2596: union sobre el tab', 'sel_qU1', 'P_U1'),
                                    ('q_D1', 'D1 SS34: sobre sus terminales', 'sel_qD1', 'P_D1'),
                                    ('q_L1', 'L1: devanado, entra por sus pads', 'sel_qL1', 'P_L1'),
                                    ('q_L2', 'L2: devanado, entra por sus pads', 'sel_qL2', 'P_L2')):
        f = ht.create(tag, 'BoundaryHeatSource', 2)
        f.label(etiqueta)
        f.selection().named(sel)
        f.set('heatSourceType', 'HeatRate')
        f.set('Pb', f'k_P*{pot}')

    if MODO == 'aire':
        aire_flujo(comp, ht, tabla)
    else:
        conveccion_formulas(comp, ht)
    operadores(comp, tabla)


def aire_flujo(comp, ht, tabla) -> None:
    """Conveccion natural calculada: flujo laminar del aire con gravedad, acoplado al calor."""
    aire = sorted(d for d, (t, _) in tabla.items() if t == 'mat_aire_cfd')
    sa = comp.selection().create('sel_aire', 'Explicit')
    sa.label('Aire')
    sa.set(aire)
    ht.feature('init1').set('Tinit', 'T_amb')
    fl = ht.create('fluido', 'FluidHeatTransferModel', 3)
    fl.label('Aire (fluido)')
    fl.selection().named('sel_aire')
    ob = ht.create('ob_T', 'OpenBoundary', 2)
    ob.label('Fronteras abiertas: entra aire a T_amb')
    ob.selection().named('sel_ext')
    for prop in ('Tustr', 'T0'):
        try:
            ob.set(prop, 'T_amb')
            break
        except Exception:
            pass
    spf = comp.physics().create('spf', 'LaminarFlow', 'geom1')
    spf.label('Flujo del aire (laminar, con gravedad)')
    spf.selection().named('sel_aire')
    pm = spf.prop('PhysicalModelProperty')
    pm.set('Compressibility', 'CompressibleMALT03')
    pm.set('IncludeGravity', '1')
    pm.set('UseReducedPressure', '1')
    # Pseudo paso de tiempo: el solver estacionario avanza como si el aire se fuera
    # asentando poco a poco; sin esto Newton no converge (se probo: 25 iteraciones y nada)
    spf.prop('AdvancedSettingProperty').set('UsePseudoTime', '1')
    try:
        ht.prop('AdvancedSettingProperty').set('UsePseudoTime', '1')
    except Exception:
        pass
    o = spf.create('ob_flujo', 'OpenBoundary', 2)
    o.label('Fronteras abiertas (sin esfuerzo normal)')
    o.selection().named('sel_ext')
    nitf = comp.multiphysics().create('nitf1', 'NonIsothermalFlow', 3)
    nitf.label('Flujo no isotermico (aire caliente sube)')
    # su seleccion es automatica: donde coinciden el flujo y el calor (el aire)


def conveccion_formulas(comp, ht) -> None:
    # Conveccion natural: la correlacion depende de hacia donde mira la cara (nz).
    # Mismas correlaciones que documenta COMSOL, escritas a mano porque las internas
    # dan NaN en este modelo.
    v = comp.variable().create('var_conv')
    v.label('Conveccion natural (correlaciones)')
    v.selection().geom('geom1', 2)
    v.selection().named('sel_ext')
    for nombre, expr, desc in (
        ('dT', 'max(T-T_amb,0.1[K])', 'Diferencia con el aire (con piso)'),
        ('beta_a', '2/(T+T_amb)', 'Expansion del aire (temperatura de pelicula)'),
        ('Ra_h', 'g_const*beta_a*dT*L_h^3/(nu_air*alpha_air)', 'Rayleigh, caras horizontales'),
        ('Ra_v', 'g_const*beta_a*dT*L_v^3/(nu_air*alpha_air)', 'Rayleigh, caras verticales'),
        ('h_arriba', 'k_air/L_h*if(Ra_h<=1e7,0.54*Ra_h^(1/4),0.15*Ra_h^(1/3))', 'Placa caliente mirando arriba'),
        ('h_abajo', 'k_air/L_h*0.27*Ra_h^(1/4)', 'Placa caliente mirando abajo'),
        ('h_vert', 'k_air/L_v*(0.825+0.387*Ra_v^(1/6)/(1+(0.492/Pr_air)^(9/16))^(8/27))^2',
         'Pared vertical (Churchill-Chu)'),
        ('h_conv', 'if(nz>0.7,h_arriba,if(nz<-0.7,h_abajo,h_vert))', 'h segun la orientacion de la cara'),
    ):
        v.set(nombre, expr, desc)
    h = ht.create('conv', 'HeatFluxBoundary', 2)
    h.label('Conveccion natural (todas las caras exteriores)')
    h.selection().named('sel_ext')
    h.set('HeatFluxType', 'ConvectiveHeatFlux')
    h.set('HeatTransferCoefficientType', 'UserDef')
    h.set('h', 'h_conv')
    h.set('Text_src', 'userdef')
    h.set('Text', 'T_amb')
    r = ht.create('rad', 'SurfaceToAmbientRadiation', 2)
    r.label('Radiacion al ambiente (apagada si eps_rad = 0)')
    r.selection().named('sel_ext')
    r.set('epsilon_rad_mat', 'userdef')
    r.set('epsilon_rad', 'eps_rad')
    r.set('Tamb_src', 'userdef')
    r.set('Tamb', 'T_amb')


def operadores(comp, tabla) -> None:
    grupos = {
        'U1': [d for d, (_, q) in tabla.items() if q in ('Cuerpo de U1', 'Tab de U1', 'Pin de U1')],
        'D1': [d for d, (_, q) in tabla.items() if q == 'Cuerpo de D1'],
        'PCB': [d for d, (_, q) in tabla.items() if q == 'Placa'],
    }
    for nombre, doms in grupos.items():
        op = comp.cpl().create(f'max_{nombre}', 'Maximum')
        op.selection().geom('geom1', 3)
        op.selection().set(doms)
    op = comp.cpl().create('max_fU1', 'Maximum')
    op.selection().geom('geom1', 2)
    op.selection().named('sel_qU1')
    op = comp.cpl().create('max_todo', 'Maximum')
    op.selection().geom('geom1', 3)
    op.selection().named('sel_todos')


def malla(comp, tabla):
    mesh = comp.mesh().create('mesh1')
    mesh.feature('size').set('hauto', '5')  # Coarse
    # Las laminas de cobre miden 35 um: se malla la geometria estirada en z y luego
    # se regresa, asi quedan elementos planos en lugar de millones de tetraedros
    ft = mesh.create('ftet1', 'FreeTet')
    ft.label('Tetraedros (geometria escalada en z)')
    ft.set('zscale', str(ESCALA_Z))
    if MODO == 'aire':
        # Con aire no hay laminas de 35 um, y un escalado fuerte deja triangulos aplastados
        # en las caras verticales que el mallador del aire no puede llenar
        ft.set('zscale', '2')
        solidos = sorted(d for d, (t, _) in tabla.items() if t != 'mat_aire_cfd')
        ft.selection().geom('geom1', 3)
        ft.selection().set(solidos)
        fa = mesh.create('ftet2', 'FreeTet')
        fa.label('Aire: tetraedros que crecen hacia afuera')
        fa.selection().geom('geom1', 3)
        fa.selection().named('sel_aire')
        sz = fa.create('size1', 'Size')
        sz.set('custom', True)
        sz.set('hmaxactive', True)
        sz.set('hmax', '15')
        sz.set('hgradactive', True)
        sz.set('hgrad', '1.5')
    mesh.run()
    return mesh


def estudios(j) -> None:
    if MODO == 'aire':
        s = j.study().create('std1')
        s.label('1 Caso 12 V -> 2 V, 0.5 A (aire simulado)')
        st = s.create('stat', 'Stationary')
        # Continuacion: se sube la potencia por etapas para que el flujo converja
        st.set('useparam', True)
        st.set('pname', ['k_P'])
        st.set('plistarr', ['0.3 1'])
        st.set('punit', [''])
        return
    s = j.study().create('std1')
    s.label('1 Caso nominal (3 A, con vias)')
    s.create('stat', 'Stationary')
    s = j.study().create('std2')
    s.label('2 Barrido de corriente')
    p = s.create('param', 'Parametric')
    p.set('pname', ['Iout'])
    p.set('plistarr', ['0.5 1 1.5 2 2.5 3'])
    p.set('punit', ['A'])
    s.create('stat', 'Stationary')
    s = j.study().create('std3')
    s.label('3 Variantes: vias, radiacion, sin L')
    p = s.create('param', 'Parametric')
    p.set('sweeptype', 'sparse')
    p.set('pname', ['con_vias', 'eps_rad', 'incluir_L'])
    p.set('plistarr', ['1 0 1 0 1', '0 0 0.9 0.9 0', '1 1 1 1 0'])
    p.set('punit', ['', '', ''])
    s.create('stat', 'Stationary')


def main() -> int:
    FIG.mkdir(exist_ok=True)
    client = mph.start(cores=4)
    try:
        model = client.load(str(ORIGEN))
        j = model.java
        comp = j.component('comp1')
        g = comp.geom('geom1')
        print('Modo:', MODO, '->', DESTINO.name)
        print('Dominios originales:', g.getNDomains())
        limpiar(j, comp)
        for nombre, expr, desc in PARAMETROS:
            j.param().set(nombre, expr, desc)
        j.param().set('k_P', '1', 'Factor de potencia (continuacion; 1 = potencia completa)')
        if MODO == 'aire':
            for nombre, valor in (('Vin', '12[V]'), ('Vout', '2[V]'), ('Iout', '0.5[A]')):
                j.param().set(nombre, valor)
        for capa, archivo in (('mF', 'mascara_cobre_F.txt'), ('mB', 'mascara_cobre_B.txt')):
            f = j.func().create(capa, 'Interpolation')
            f.label(f'Mascara de cobre {capa[1]}.Cu (1 = cobre)')
            f.set('source', 'file')
            f.set('filename', str(GEO / archivo))
            f.set('struct', 'grid')
            f.set('nargs', '2')
            f.importData()
            f.set('fununit', '1')
            f.set('extrap', 'const')
        preparar_geometria(g)
        print('Dominios del modelo:', g.getNDomains())
        tabla = clasificar(g)
        conteo = {}
        for t, _ in tabla.values():
            conteo[t] = conteo.get(t, 0) + 1
        print('Dominios por material:', conteo)
        with (FEM / 'geometria' / f'dominios_materiales_STEP_{MODO}.csv').open('w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['dominio', 'material', 'que es'])
            etiquetas = {m[0]: m[1] for m in MATERIALES}
            for d in sorted(tabla):
                w.writerow([d, etiquetas[tabla[d][0]], tabla[d][1]])
        materiales(comp, tabla)
        fisica(comp, tabla, g.getNDomains())
        estudios(j)
        model.save(str(DESTINO))
        print('Guardado (sin malla):', DESTINO)
        mesh = malla(comp, tabla)
        print('Elementos de malla:', mesh.getNumElem())
        model.save(str(DESTINO))
        print('Guardado (con malla):', DESTINO)
        return 0
    finally:
        client.disconnect()


if __name__ == '__main__':
    raise SystemExit(main())
