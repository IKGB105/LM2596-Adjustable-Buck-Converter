"""Resuelve el modelo con aire simulado (MODO=aire de modelo_termico_step.py).

Caso 12 V -> 2 V, 0.5 A. La conveccion natural se CALCULA (flujo laminar del aire
con gravedad, acoplado al calor), en lugar de usar correlaciones.

Salidas:
  - resultados/caso_12V_2V_05A_aire.csv
  - figuras/aire_temperatura_placa.png      (temperatura en la placa y componentes)
  - figuras/aire_corte_temperatura.png      (corte vertical: aire caliente subiendo)
  - figuras/aire_corte_velocidad.png        (corte vertical: velocidad del aire)
"""
from __future__ import annotations

import csv
from pathlib import Path

import mph
import numpy as np

FEM = Path(__file__).resolve().parent.parent
MODELO = FEM / 'BUCK_LM2596_termico_12V_2V_05A_aire.mph'
TEMPS = [('T union U1 (cara tab-cuerpo)', 'max_fU1(T)'), ('T max U1', 'max_U1(T)'),
         ('T max D1', 'max_D1(T)'), ('T max placa', 'max_PCB(T)')]
Y_CORTE = -25.5  # mm: plano vertical que pasa por el tab de U1


def exportar(j, tag, archivo):
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


def grupo(j, tag, titulo, dset):
    if tag in [str(t) for t in j.result().tags()]:
        j.result().remove(tag)
    pg = j.result().create(tag, 'PlotGroup3D')
    pg.set('data', dset)
    pg.set('titletype', 'manual')
    pg.set('title', titulo)
    pg.label(titulo[:60])
    return pg


def main() -> int:
    client = mph.start(cores=4)
    try:
        model = client.load(str(MODELO))
        j = model.java
        comp = j.component('comp1')
        # La seleccion del aire puede llegar vacia al recargar el modelo: se vuelve a fijar
        # ubicando el dominio por tamano (la caja de aire mide mas de 100 mm)
        g = comp.geom('geom1')
        mf = g.measureFinal()
        mf.selection().geom(3)
        aire = []
        for d_ in range(1, g.getNDomains() + 1):
            mf.selection().set([d_])
            bb = list(mf.getBoundingBox())
            if bb[1] - bb[0] > 90:
                aire.append(d_)
        comp.selection('sel_aire').set(aire)
        for nombre, ents in (('sel_aire', comp.selection('sel_aire').entities()),
                             ('flujo spf', comp.physics('spf').selection().entities()),
                             ('fluido (calor)', comp.physics('ht').feature('fluido').selection().entities()),
                             ('acoplamiento nitf', comp.multiphysics('nitf1').selection().entities())):
            ents = [int(x) for x in ents]
            print(f'  {nombre}: dominios {ents}', flush=True)
            assert ents, f'{nombre} quedo sin dominios: no se simularia el aire'
        # Las selecciones de la malla tambien llegan mal: se fijan y se vuelve a mallar
        solidos_m = [d_ for d_ in range(1, g.getNDomains() + 1) if d_ not in aire]
        mesh = comp.mesh('mesh1')
        f1 = mesh.feature('ftet1')
        f1.selection().geom('geom1', 3)
        f1.selection().set(solidos_m)
        f2 = mesh.feature('ftet2')
        f2.selection().geom('geom1', 3)
        f2.selection().set(aire)
        print(f'  malla: solidos {len(list(f1.selection().entities()))} dominios, '
              f'aire {list(f2.selection().entities())}', flush=True)
        mesh.run()
        print(f'  malla: {mesh.getNumElem()} elementos, calidad min {mesh.stat().getMinQuality():.3f}', flush=True)
        model.save(str(MODELO))
        # Operador de maximo sobre el aire: se crea ANTES de resolver para que la solucion lo conozca
        if 'max_aire' not in [str(t) for t in comp.cpl().tags()]:
            op = comp.cpl().create('max_aire', 'Maximum')  # la velocidad solo existe en el aire
            op.selection().geom('geom1', 3)
            op.selection().named('sel_aire')
        etiqueta = str(j.study('std1').label())
        # Secuencia de solver explicita para subir el limite de iteraciones de Newton
        # (con 25, el valor por defecto, no convergio)
        for t in [str(x) for x in j.sol().tags()]:
            j.sol().remove(t)
        sol = j.sol().create('sol1')
        sol.study('std1')
        sol.attach('std1')
        sol.createAutoSequence('std1')
        s1 = sol.feature('s1')
        for t in [str(x) for x in s1.feature().tags()]:
            f = s1.feature(t)
            tipo = str(f.getType())
            for prop, val in (('maxiter', '200'), ('maxsegiter', '200')):
                try:
                    f.set(prop, val)
                    print(f'  {t} ({tipo}): {prop} = {val}', flush=True)
                except Exception:
                    pass
        print('Resolviendo:', etiqueta, flush=True)
        sol.runAll()
        d = next(x for x in model / 'datasets' if x.tag() == 'dset1')
        n_pasos = 2  # k_P = 0.3 1 -> el ultimo es la potencia completa

        def ev(expr, unidad=None):
            return float(np.real(np.ravel(model.evaluate(expr, unidad, dataset=d, inner=[n_pasos]))[0]))

        filas = [('Vin', ev('Vin', 'V'), 'V'), ('Vout', ev('Vout', 'V'), 'V'), ('Iout', ev('Iout', 'A'), 'A'),
                 ('Factor de potencia k_P (1 = completa)', ev('k_P'), '')]
        filas += [(p, ev(p, 'W'), 'W') for p in ('P_U1', 'P_D1', 'P_L1', 'P_L2', 'P_tot')]
        filas += [(n, ev(e, 'degC'), 'degC') for n, e in TEMPS]
        filas += [('Velocidad maxima del aire', ev('max_aire(spf.U)', 'mm/s'), 'mm/s'),
                  ('T max del aire', ev('max_aire(T)', 'degC'), 'degC')]
        (FEM / 'resultados').mkdir(exist_ok=True)
        with (FEM / 'resultados' / 'caso_12V_2V_05A_aire.csv').open('w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['cantidad', 'valor', 'unidad'])
            for n, v, u in filas:
                w.writerow([n, f'{v:.4g}', u])
                print(f'  {n}: {v:.4g} {u}', flush=True)

        # Dataset solo con los solidos (para ver la placa sin la caja de aire encima)
        solidos = [d_ for d_ in range(1, comp.geom('geom1').getNDomains() + 1)
                   if d_ not in [int(x) for x in comp.selection('sel_aire').entities()]]
        if 'dset_sol' in [str(t) for t in j.result().dataset().tags()]:
            j.result().dataset().remove('dset_sol')
        ds = j.result().dataset().create('dset_sol', 'Solution')
        ds.set('solution', 'sol1')
        ds.selection().geom('geom1', 3)
        ds.selection().set(solidos)

        pg = grupo(j, 'pg_aire_T', 'Temperatura [degC] con aire simulado: 12 V -> 2 V, 0.5 A, T amb 25 degC', 'dset_sol')
        pg.setIndex('looplevel', str(n_pasos), 0)
        s = pg.create('surf1', 'Surface')
        s.set('expr', 'T')
        s.set('unit', 'degC')
        s.set('colortable', 'HeatCameraLight')
        exportar(j, 'pg_aire_T', FEM / 'figuras' / 'aire_temperatura_placa.png')

        for tag, expr, unidad, titulo, archivo, tabla_color in (
            ('pg_aire_cT', 'T', 'degC', 'Corte vertical por U1: temperatura del aire [degC] (el aire caliente sube)',
             'aire_corte_temperatura.png', 'HeatCameraLight'),
            ('pg_aire_cU', 'spf.U', 'mm/s', 'Corte vertical por U1: velocidad del aire [mm/s] (conveccion natural)',
             'aire_corte_velocidad.png', 'Rainbow'),
        ):
            pg = grupo(j, tag, titulo, 'dset1')
            pg.setIndex('looplevel', str(n_pasos), 0)
            sl = pg.create('slc1', 'Slice')
            sl.set('expr', expr)
            sl.set('unit', unidad)
            sl.set('quickplane', 'zx')
            sl.set('quickynumber', '1')
            sl.set('quickymethod', 'coord')
            sl.set('quicky', f'{Y_CORTE:g}')
            sl.set('colortable', tabla_color)
            sp = pg.create('surf_sol', 'Surface')  # silueta de la placa
            sp.set('data', 'dset_sol')
            sp.set('expr', 'T')
            sp.set('unit', 'degC')
            sp.set('coloring', 'uniform')
            sp.set('color', 'gray')
            exportar(j, tag, FEM / 'figuras' / archivo)

        model.save(str(MODELO))
        print('Guardado:', MODELO)
        return 0
    finally:
        client.disconnect()


if __name__ == '__main__':
    raise SystemExit(main())
