"""Simulacion termica (COMSOL, Heat Transfer in Solids) de la placa BUCK LM2596.

Construye fem/BUCK_LM2596_termico.mph con geometria primitiva (sin importar STEP):
  - Placa FR4 60 x 30 x 1.6 mm.
  - U1 LM2596S (TO-263): cuerpo EMC + tab de cobre (solido separado).
  - D1 SS34 (SMC), L1 33 uH y L2 20 uH (12 x 12 x 4.5 mm), C1-C3 como cilindros.
  - 16 vias (cilindros de cobre efectivo) que se pueden apagar con `con_vias`.
  - Cobre de F.Cu y B.Cu como capas delgadas de 35 um ("Thin Layer"), con la
    forma real de pads, pistas y plano GND (mascaras exportadas de KiCad).
Las medidas y posiciones vienen del STEP que exporta KiCad (geometria/).

Convencion de coordenadas: X = x_KiCad - 20, Y = 50 - y_KiCad (mm), Z hacia arriba.
Igual que en FuenteLM317_9V: sin dominio de aire, conveccion natural con la
correlacion de cada tipo de cara y potencias como "Heat rate".
"""
from __future__ import annotations

from pathlib import Path

import mph
import numpy as np

FEM = Path(__file__).resolve().parent.parent
GEO = FEM / 'geometria'
FIG = FEM / 'figuras'
T_AMB = 25.0

# Vias (coordenadas KiCad, mm) -> X = x - 20, Y = 50 - y
VIAS_KICAD = [(35.95, 28.8), (31.65, 22.4), (35.9, 22.4), (31.6, 25.3), (35.85, 25.3), (40.2, 28.8),
              (39.2, 32.2), (48.7, 30.0), (40.1, 25.3), (35.8, 32.3), (31.7, 28.8), (26.6, 42.8),
              (75.6, 47.8), (24.7, 47.8), (40.15, 22.4), (25.4, 35.7)]


def parametros(p) -> None:
    datos = [
        # Punto de operacion
        ('Vin', '35[V]', 'Voltaje de entrada (esquematico: 35 V -> 1.2 V)'),
        ('Vout', '1.23[V]', 'Voltaje de salida (minimo del LM2596-ADJ)'),
        ('Iout', '3[A]', 'Corriente de carga (peor caso: maximo del LM2596)'),
        ('fsw', '150[kHz]', 'Frecuencia de conmutacion del LM2596 (datasheet)'),
        # Supuestos de perdidas
        ('Vsat', '1.16[V]', 'Caida del switch interno del LM2596 a 3 A (datasheet, tipico)'),
        ('Iq', '5[mA]', 'Corriente de reposo del LM2596 (datasheet, tipico)'),
        ('t_sw', '100[ns]', 'Tiempo de conmutacion total subida+bajada (SUPUESTO)'),
        ('Vf', '0.5[V]', 'Caida del diodo SS34 a 3 A (datasheet, maximo)'),
        ('DCR1', '40[mohm]', 'Resistencia del devanado de L1 33 uH 12x12 (SUPUESTO)'),
        ('DCR2', '25[mohm]', 'Resistencia del devanado de L2 20 uH 12x12 (SUPUESTO)'),
        ('incluir_L', '1', '1 = incluir perdidas de L1 y L2, 0 = ignorarlas'),
        ('con_vias', '1', '1 = vias de cobre, 0 = sin vias (solo FR4)'),
        ('eps_rad', '0', 'Emisividad para radiacion (0 = solo conveccion, como en la LM317)'),
        ('T_amb', f'{T_AMB}[degC]', 'Temperatura ambiente'),
        ('theta_jc', '2[K/W]', 'Resistencia union-tab del LM2596 TO-263 (datasheet)'),
        # Perdidas calculadas
        ('Dcic', '(Vout+Vf)/(Vin-Vsat+Vf)', 'Ciclo de trabajo'),
        ('P_U1', 'Vsat*Iout*Dcic+0.5*Vin*Iout*t_sw*fsw+Vin*Iq', 'Perdidas del LM2596'),
        ('P_D1', 'Vf*Iout*(1-Dcic)', 'Perdidas del diodo SS34'),
        ('P_L1', 'incluir_L*Iout^2*DCR1', 'Perdidas de L1'),
        ('P_L2', 'incluir_L*Iout^2*DCR2', 'Perdidas de L2'),
        ('P_tot', 'P_U1+P_D1+P_L1+P_L2', 'Perdidas totales'),
        ('eficiencia', 'Vout*Iout/(Vout*Iout+P_tot)', 'Eficiencia estimada'),
        # Materiales
        ('k_cu', '400[W/(m*K)]', 'Conductividad del cobre'),
        ('t_cu', '35[um]', 'Espesor del cobre (1 oz)'),
        ('k_via', '122[W/(m*K)]', 'k efectiva de una via: barril de 25 um en taladro de 0.3 mm'),
        ('k_fr4', '0.3[W/(m*K)]', 'Conductividad del FR4'),
        # Aire (para la correlacion de las caras horizontales hacia arriba)
        ('k_air', '0.028[W/(m*K)]', 'Conductividad del aire (~50 degC)'),
        ('nu_air', '1.8e-5[m^2/s]', 'Viscosidad cinematica del aire (~50 degC)'),
        ('alpha_air', '2.55e-5[m^2/s]', 'Difusividad termica del aire (~50 degC)'),
        ('L_up', '10[mm]', 'Longitud caracteristica de las caras de arriba (area/perimetro de la placa)'),
    ]
    for nombre, expr, desc in datos:
        p.set(nombre, expr, desc)


def bloque(geom, tag, etiqueta, x, y, z):
    b = geom.create(tag, 'Block')
    b.label(etiqueta)
    b.set('pos', [f'{x[0]:g}', f'{y[0]:g}', f'{z[0]:g}'])
    b.set('size', [f'{x[1] - x[0]:g}', f'{y[1] - y[0]:g}', f'{z[1] - z[0]:g}'])
    b.set('selresult', True)


def cilindro(geom, tag, etiqueta, cx, cy, r, z0, h):
    c = geom.create(tag, 'Cylinder')
    c.label(etiqueta)
    c.set('pos', [f'{cx:g}', f'{cy:g}', f'{z0:g}'])
    c.set('r', f'{r:g}')
    c.set('h', f'{h:g}')
    c.set('selresult', True)


def caja_sel(geom, tag, dim, x, y, z, condicion='inside'):
    b = geom.create(tag, 'BoxSelection')
    b.set('entitydim', str(dim))
    b.set('condition', condicion)
    for k, v in zip(('xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax'), (*x, *y, *z)):
        b.set(k, f'{v:g}')


def geometria(comp):
    geom = comp.geom().create('geom1', 3)
    geom.lengthUnit('mm')
    t = 1.6  # superficie de la placa
    bloque(geom, 'pcb', 'Placa FR4', (0, 60), (0, 30), (0, t))
    bloque(geom, 'u1', 'U1 cuerpo (EMC)', (10.8, 20.8), (18.07, 27.33), (t, t + 4.42))
    bloque(geom, 'tab', 'U1 tab (cobre)', (10.8, 20.8), (20.77, 28.32), (t, t + 1.27))
    bloque(geom, 'd1', 'D1 SS34 (SMC)', (8.1, 16.1), (4.25, 10.15), (t, t + 2.31))
    bloque(geom, 'l1', 'L1 33 uH', (22.5, 34.5), (1.9, 13.9), (t, t + 4.5))
    bloque(geom, 'l2', 'L2 20 uH', (38.75, 50.75), (2.1, 14.1), (t, t + 4.5))
    cilindro(geom, 'c1', 'C1 100 uF', 5.4, 17.7, 4.0, t, 8.0)
    cilindro(geom, 'c2', 'C2 220 uF', 38.0, 20.9, 5.0, t, 10.0)
    cilindro(geom, 'c3', 'C3 100 uF', 48.15, 20.7, 4.0, t, 8.0)
    for i, (xk, yk) in enumerate(VIAS_KICAD, start=1):
        cilindro(geom, f'via{i}', f'Via {i}', xk - 20, 50 - yk, 0.15, 0, t)
    u = geom.create('uni_vias', 'Union')
    u.label('Vias (grupo)')
    u.selection('input').set([f'via{i}' for i in range(1, len(VIAS_KICAD) + 1)])
    u.set('intbnd', True)
    u.set('selresult', True)
    geom.feature('fin').set('action', 'union')

    # Caras por posicion (mm)
    caja_sel(geom, 'sel_cara_sup', 2, (-1, 61), (-1, 31), (t - 0.01, t + 0.01))
    caja_sel(geom, 'sel_cara_inf', 2, (-1, 61), (-1, 31), (-0.01, 0.01))
    # Fuentes: caras que cruzan una caja delgada en el centro de la pieza (ahi solo
    # hay caras horizontales), mas robusto que exigir la cara completa dentro
    caja_sel(geom, 'sel_fuente_u1', 2, (12, 19), (22, 26), (t + 1.26, t + 1.28), 'intersects')
    caja_sel(geom, 'sel_fuente_d1', 2, (9.5, 14.5), (5.5, 9), (t - 0.01, t + 0.01), 'intersects')
    # Caras horizontales hacia arriba (techo de cada pieza y placa expuesta)
    alturas = [t, t + 1.27, t + 2.31, t + 4.42, t + 4.5, t + 8.0, t + 10.0]
    for i, zz in enumerate(alturas):
        caja_sel(geom, f'sel_h{i}', 2, (-1, 61), (-1, 31), (zz - 0.005, zz + 0.005))
    un = geom.create('sel_horiz', 'UnionSelection')
    un.set('entitydim', '2')
    un.set('input', [f'sel_h{i}' for i in range(len(alturas))])
    geom.run()
    return geom


def materiales(comp):
    datos = [
        ('mat_fr4', 'FR4', 'geom1_pcb_dom', 'k_fr4', '1900', '1369'),
        ('mat_emc', 'EMC (cuerpo de U1 y D1)', None, '0.7', '1900', '900'),
        ('mat_fer', 'Ferrita (L1, L2)', None, '4', '4800', '750'),
        ('mat_abs', 'Relleno generico (capacitores)', None, '0.2', '1050', '1400'),
        ('mat_cu', 'Cobre (tab de U1)', 'geom1_tab_dom', 'k_cu', '8960', '385'),
        ('mat_via', 'Vias (cobre efectivo)', 'geom1_uni_vias_dom', 'if(con_vias,k_via,k_fr4)', '8960', '385'),
    ]
    sel_extra = {
        'mat_emc': ['geom1_u1_dom', 'geom1_d1_dom'],
        'mat_fer': ['geom1_l1_dom', 'geom1_l2_dom'],
        'mat_abs': ['geom1_c1_dom', 'geom1_c2_dom', 'geom1_c3_dom'],
    }
    for tag, etiqueta, sel, k, rho, cp in datos:
        m = comp.material().create(tag, 'Common')
        m.label(etiqueta)
        if sel:
            m.selection().named(sel)
        else:
            u = comp.selection().create(f'sel_{tag}', 'Union')
            u.set('input', sel_extra[tag])
            m.selection().named(f'sel_{tag}')
        pg = m.propertyGroup('def')
        pg.set('thermalconductivity', [k])
        pg.set('density', [rho])
        pg.set('heatcapacity', [cp])


def operadores(comp):
    for tag, sel, dim in (('max_U1', ['geom1_u1_dom', 'geom1_tab_dom'], 3),
                          ('max_D1', ['geom1_d1_dom'], 3), ('max_L1', ['geom1_l1_dom'], 3),
                          ('max_L2', ['geom1_l2_dom'], 3), ('max_C1', ['geom1_c1_dom'], 3),
                          ('max_C2', ['geom1_c2_dom'], 3), ('max_PCB', ['geom1_pcb_dom'], 3),
                          ('max_fU1', ['geom1_sel_fuente_u1'], 2)):
        op = comp.cpl().create(tag, 'Maximum')
        op.selection().geom('geom1', dim)
        if len(sel) == 1:
            op.selection().named(sel[0])
        else:
            u = comp.selection().create(f'sel_{tag}', 'Union')
            u.set('input', sel)
            op.selection().named(f'sel_{tag}')


def fisica(comp, j):
    for capa, archivo in (('mF', 'mascara_cobre_F.txt'), ('mB', 'mascara_cobre_B.txt')):
        f = j.func().create(capa, 'Interpolation')
        f.label(f'Mascara de cobre {capa[1]}.Cu')
        f.set('source', 'file')
        f.set('filename', str(GEO / archivo))
        f.set('struct', 'grid')
        f.set('nargs', '2')
        f.importData()
        f.set('argunit', ['mm', 'mm'])
        f.set('fununit', '1')
        f.set('extrap', 'const')

    ht = comp.physics().create('ht', 'HeatTransfer', 'geom1')
    # Arranque 20 K arriba del ambiente: con T = T_amb exacto la correlacion de
    # conveccion natural tiene derivada infinita en (T - T_amb)^(1/4)
    ht.feature('init1').set('Tinit', 'T_amb+20[K]')

    # Cobre de cada cara como capa delgada conductiva de 35 um
    for tag, etiqueta, sel, mascara in (('tl_F', 'Cobre F.Cu (pads y pistas)', 'geom1_sel_cara_sup', 'mF'),
                                        ('tl_B', 'Cobre B.Cu (plano GND)', 'geom1_sel_cara_inf', 'mB')):
        tl = ht.create(tag, 'ThinLayer', 2)
        tl.label(etiqueta)
        tl.selection().named(sel)
        tl.set('LayerType', 'Conductive')
        tl.set('ds', 't_cu')
        tl.set('ks_mat', 'userdef')
        k = f'k_fr4+(k_cu-k_fr4)*{mascara}(x,y)'
        tl.set('ks', [k, '0', '0', '0', k, '0', '0', '0', k])
        for prop, val in (('rhos_mat', 'userdef'), ('rhos', '8960'), ('Cp_s_mat', 'userdef'), ('Cp_s', '385')):
            tl.set(prop, val)

    # Fuentes de calor (modo Heat rate: potencia total en W)
    f = ht.create('bhs_U1', 'BoundaryHeatSource', 2)
    f.label('U1: union sobre el tab')
    f.selection().named('geom1_sel_fuente_u1')
    f.set('heatSourceType', 'HeatRate')
    f.set('Pb', 'P_U1')
    f = ht.create('bhs_D1', 'BoundaryHeatSource', 2)
    f.label('D1: union sobre sus terminales')
    f.selection().named('geom1_sel_fuente_d1')
    f.set('heatSourceType', 'HeatRate')
    f.set('Pb', 'P_D1')
    for tag, sel, pot in (('hs_L1', 'geom1_l1_dom', 'P_L1'), ('hs_L2', 'geom1_l2_dom', 'P_L2')):
        f = ht.create(tag, 'HeatSource', 3)
        f.label(f'{tag[3:]}: perdidas del devanado')
        f.selection().named(sel)
        f.set('heatSourceType', 'HeatRate')
        f.set('P0', pot)

    # Conveccion natural por tipo de cara (como en la LM317)
    ext = comp.selection().create('sel_exterior', 'Adjacent')
    ext.label('Todas las caras exteriores')
    ext.set('input', ['geom1_pcb_dom', 'geom1_u1_dom', 'geom1_tab_dom', 'geom1_d1_dom', 'geom1_l1_dom',
                      'geom1_l2_dom', 'geom1_c1_dom', 'geom1_c2_dom', 'geom1_c3_dom', 'geom1_uni_vias_dom'])
    arriba = comp.selection().create('sel_arriba', 'Intersection')
    arriba.set('entitydim', '2')
    arriba.set('input', ['sel_exterior', 'geom1_sel_horiz'])
    abajo = comp.selection().create('sel_abajo', 'Intersection')
    abajo.set('entitydim', '2')
    abajo.set('input', ['sel_exterior', 'geom1_sel_cara_inf'])
    lados = comp.selection().create('sel_lados', 'Difference')
    lados.set('entitydim', '2')
    lados.set('add', ['sel_exterior'])
    lados.set('subtract', ['sel_arriba', 'sel_abajo'])
    for tag, etiqueta, sel, tipo, longitud in (
        ('hf_abajo', 'Conveccion: cara inferior de la placa', 'sel_abajo', 'HorizontalPlateDown', ('Lpl', '10[mm]')),
        ('hf_lados', 'Conveccion: caras verticales', 'sel_lados', 'VerticalWall', ('Lwall', '5[mm]')),
    ):
        h = ht.create(tag, 'HeatFluxBoundary', 2)
        h.label(etiqueta)
        h.selection().named(sel)
        h.set('HeatFluxType', 'ConvectiveHeatFlux')
        h.set('HeatTransferCoefficientType', 'ExtNaturalConvection')
        h.set('ExtNaturalConvectionType', tipo)
        h.set(*longitud)
        h.set('Text_src', 'userdef')
        h.set('Text', 'T_amb')
        h.set('pA_src', 'userdef')
        h.set('pA', '1[atm]')
    # Caras hacia arriba: la correlacion interna "Horizontal plate, upside" da NaN en
    # este modelo, asi que se escribe la misma correlacion de forma explicita:
    # Nu = 0.54 Ra^(1/4) (Ra <= 1e7), 0.15 Ra^(1/3) (Ra > 1e7), aire a T de pelicula
    v = comp.variable().create('var_conv')
    v.label('Conveccion natural en caras horizontales hacia arriba')
    v.selection().geom('geom1', 2)
    v.selection().named('sel_arriba')
    v.set('dT_up', 'max(T-T_amb,0.1[K])', 'Diferencia de temperatura con el aire (con piso)')
    v.set('Ra_up', 'g_const*(2/(T+T_amb))*dT_up*L_up^3/(nu_air*alpha_air)', 'Numero de Rayleigh')
    v.set('h_up', 'k_air/L_up*if(Ra_up<=1e7,0.54*Ra_up^(1/4),0.15*Ra_up^(1/3))', 'Coeficiente de conveccion')
    h = ht.create('hf_arriba', 'HeatFluxBoundary', 2)
    h.label('Conveccion: caras horizontales hacia arriba (correlacion explicita)')
    h.selection().named('sel_arriba')
    h.set('HeatFluxType', 'ConvectiveHeatFlux')
    h.set('HeatTransferCoefficientType', 'UserDef')
    h.set('h', 'h_up')
    h.set('Text_src', 'userdef')
    h.set('Text', 'T_amb')
    r = ht.create('rad', 'SurfaceToAmbientRadiation', 2)
    r.label('Radiacion al ambiente (apagada si eps_rad = 0)')
    r.selection().named('sel_exterior')
    r.set('epsilon_rad_mat', 'userdef')
    r.set('epsilon_rad', 'eps_rad')
    r.set('Tamb_src', 'userdef')
    r.set('Tamb', 'T_amb')


def malla(comp):
    mesh = comp.mesh().create('mesh1')
    mesh.autoMeshSize(4)  # Normal
    mesh.run()
    return mesh


def estudios(j):
    s1 = j.study().create('std1')
    s1.label('1 Caso nominal (3 A, con vias)')
    s1.create('stat', 'Stationary')

    s2 = j.study().create('std2')
    s2.label('2 Barrido de corriente')
    p = s2.create('param', 'Parametric')
    p.set('pname', ['Iout'])
    p.set('plistarr', ['0.5 1 1.5 2 2.5 3'])
    p.set('punit', ['A'])
    s2.create('stat', 'Stationary')

    s3 = j.study().create('std3')
    s3.label('3 Variantes: vias, radiacion, sin L')
    p = s3.create('param', 'Parametric')
    p.set('sweeptype', 'sparse')
    p.set('pname', ['con_vias', 'eps_rad', 'incluir_L'])
    p.set('plistarr', ['1 0 1 0 1', '0 0 0.9 0.9 0', '1 1 1 1 0'])
    p.set('punit', ['', '', ''])
    s3.create('stat', 'Stationary')


VARIANTES = ['Con vias (base)', 'Sin vias', 'Con vias + radiacion', 'Sin vias + radiacion', 'Con vias, sin perdidas en L']
EXPR = ['max_fU1(T)', 'max_U1(T)', 'max_D1(T)', 'max_L1(T)', 'max_L2(T)', 'max_C1(T)', 'max_C2(T)', 'max_PCB(T)']


def graficas(model, j):
    pg = j.result().create('pg_T', 'PlotGroup3D')
    pg.label('Temperatura (caso nominal)')
    pg.set('data', 'dset1')
    pg.set('titletype', 'manual')
    pg.set('title', 'Temperatura [degC]: Vin 35 V, Vout 1.23 V, Iout 3 A, con vias, solo conveccion')
    s = pg.create('surf1', 'Surface')
    s.set('expr', 'T')
    s.set('unit', 'degC')
    s.set('colortable', 'HeatCameraLight')

    pc = j.result().create('pg_corte', 'PlotGroup3D')
    pc.label('Temperatura en la cara superior de la placa')
    pc.set('data', 'dset1')
    pc.set('titletype', 'manual')
    pc.set('title', 'Temperatura [degC] en el plano de la placa (z = 1.6 mm): como se reparte el calor por el cobre')
    sl = pc.create('slc1', 'Slice')
    sl.set('expr', 'T')
    sl.set('unit', 'degC')
    sl.set('quickplane', 'xy')
    sl.set('quickzmethod', 'coord')
    sl.set('quickz', '1.55')
    sl.set('colortable', 'HeatCameraLight')


def exportar(model, tag, archivo, vista_superior=False):
    j = model.java
    e = j.result().export().create('img_tmp', 'Image')
    e.set('plotgroup', tag)
    e.set('filename', str(archivo))
    e.set('size', 'manualweb')
    e.set('width', '1400')
    e.set('height', '800')
    e.set('options3d', 'on')
    for p in ('title3d', 'legend3d'):
        e.set(p, 'on')
    e.set('background', 'color')
    try:
        e.run()
    finally:
        j.result().export().remove('img_tmp')


def main() -> int:
    FIG.mkdir(parents=True, exist_ok=True)
    client = mph.start(cores=4)
    try:
        model = client.create('BUCK_LM2596_termico')
        j = model.java
        parametros(j.param())
        comp = j.component().create('comp1', True)
        geometria(comp)
        materiales(comp)
        operadores(comp)
        fisica(comp, j)
        for sel in ('geom1_sel_fuente_u1', 'geom1_sel_fuente_d1'):
            n = len(list(comp.selection(sel).entities()))
            print(f'{sel}: {n} cara(s)')
            assert n >= 1, f'La seleccion {sel} quedo vacia'
        mesh = malla(comp)
        print('Elementos de malla:', mesh.getNumElem())
        estudios(j)
        guardar = FEM / 'BUCK_LM2596_termico.mph'
        model.save(str(guardar))  # guardar antes de resolver por si algo falla

        print('Resolviendo caso nominal...')
        model.solve('1 Caso nominal (3 A, con vias)')
        pot = {n: float(np.ravel(model.evaluate(n, 'W'))[0]) for n in ('P_U1', 'P_D1', 'P_L1', 'P_L2', 'P_tot')}
        efi = float(np.ravel(model.evaluate('eficiencia'))[0])
        tn = [float(np.ravel(model.evaluate(e, 'degC'))[0]) for e in EXPR]
        print('Potencias [W]:', {k: round(v, 3) for k, v in pot.items()}, 'eficiencia', round(efi, 3))
        print('T max [degC] fuenteU1, U1, D1, L1, L2, C1, C2, PCB:', [round(v, 1) for v in tn])
        graficas(model, j)
        exportar(model, 'pg_T', FIG / 'temperatura_3D.png')
        exportar(model, 'pg_corte', FIG / 'temperatura_plano_placa.png')
        model.save(str(guardar))

        print('Resolviendo barrido de corriente...')
        model.solve('2 Barrido de corriente')
        print('Resolviendo variantes...')
        model.solve('3 Variantes: vias, radiacion, sin L')
        model.save(str(guardar))
        print('Guardado:', guardar)
        return 0
    finally:
        client.disconnect()


if __name__ == '__main__':
    raise SystemExit(main())
