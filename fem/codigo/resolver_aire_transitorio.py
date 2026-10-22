"""Aire simulado por TRANSITORIO hasta estado estable (12 V -> 2 V, 0.5 A).

El estacionario de conveccion natural acoplada no convergio (Newton, 200 iteraciones,
con pseudo paso de tiempo). Aqui se simula el calentamiento en el tiempo hasta que la
temperatura se estabiliza.

Truco valido para llegar rapido al estado estable: la densidad de los SOLIDOS se
multiplica por esc_C = 0.01 (capacidad termica 100 veces menor). La temperatura de
equilibrio no depende de la capacidad termica, solo cuanto tarda en llegar. El aire
conserva su densidad real (es la que produce la flotacion).

Al final se verifica que T ya no cambie (ultimos dos tiempos).
"""
from __future__ import annotations

import csv
from pathlib import Path

import mph
import numpy as np

FEM = Path(__file__).resolve().parent.parent
MODELO = FEM / 'BUCK_LM2596_termico_12V_2V_05A_aire.mph'
T_FIN = 300  # s (con capacidad 100x menor, la placa se estabiliza en unos ~15 s)
TEMPS = [('T union U1 (cara tab-cuerpo)', 'max_fU1(T)'), ('T max U1', 'max_U1(T)'),
         ('T max D1', 'max_D1(T)'), ('T max placa', 'max_PCB(T)')]


def main() -> int:
    client = mph.start(cores=4)
    try:
        model = client.load(str(MODELO))
        j = model.java
        comp = j.component('comp1')
        g = comp.geom('geom1')

        # Verificaciones (las selecciones ya se corrigieron en resolver_aire.py y se guardaron)
        aire = [int(x) for x in comp.selection('sel_aire').entities()]
        for nombre, ents in (('flujo spf', comp.physics('spf').selection().entities()),
                             ('acoplamiento nitf', comp.multiphysics('nitf1').selection().entities())):
            print(f'  {nombre}: dominios {[int(x) for x in ents]}', flush=True)
            assert list(ents), f'{nombre} sin dominios'
        assert aire, 'sel_aire vacia'

        # Capacidad termica de los solidos reducida (no la del aire)
        j.param().set('esc_C', '0.01', 'Factor de capacidad termica de solidos (solo acelera el transitorio)')
        j.param().set('k_P', '1')
        for t in [str(x) for x in comp.material().tags()]:
            if t == 'mat_aire_cfd':
                continue
            pg = comp.material(t).propertyGroup('def')
            rho = str(pg.getStringArray('density')[0])
            if 'esc_C' not in rho:
                pg.set('density', [f'({rho})*esc_C'])
        if 'max_aire' not in [str(t) for t in comp.cpl().tags()]:
            op = comp.cpl().create('max_aire', 'Maximum')
            op.selection().geom('geom1', 3)
            op.selection().set(aire)

        if 'std_t' in [str(t) for t in j.study().tags()]:
            j.study().remove('std_t')
        s = j.study().create('std_t')
        s.label('2 Transitorio hasta estado estable (aire simulado)')
        td = s.create('time', 'Transient')
        td.set('tlist', f'range(0,10,{T_FIN})')
        etiqueta = str(s.label())
        print('Resolviendo:', etiqueta, flush=True)
        antes = {d.tag() for d in model / 'datasets'}
        model.solve(etiqueta)
        d = [x for x in model / 'datasets' if x.tag() not in antes][-1]
        model.save(str(MODELO))  # guardar ANTES de evaluar: si algo falla, la solucion no se pierde
        print('Solucion guardada:', MODELO, flush=True)

        def ev(expr, unidad=None, paso='last'):
            return float(np.real(np.ravel(model.evaluate(expr, unidad, dataset=d, inner=paso))[-1]))

        filas = [('Vin', ev('Vin', 'V'), 'V'), ('Vout', ev('Vout', 'V'), 'V'), ('Iout', ev('Iout', 'A'), 'A')]
        filas += [(p, ev(p, 'W'), 'W') for p in ('P_U1', 'P_D1', 'P_L1', 'P_L2', 'P_tot')]
        filas += [(n, ev(e, 'degC'), 'degC') for n, e in TEMPS]
        filas += [('Velocidad maxima del aire', ev('max_aire(spf.U)', 'mm/s'), 'mm/s'),
                  ('T max del aire', ev('max_aire(T)', 'degC'), 'degC')]
        # Estabilidad: cambio de la T de U1 entre los dos ultimos tiempos guardados
        tiempos = np.ravel(model.evaluate('t', 's', dataset=d))  # sin "inner": todos los tiempos
        t_u1 = np.ravel(model.evaluate('max_fU1(T)', 'degC', dataset=d))
        if len(t_u1) >= 2:
            filas.append(('Cambio de T U1 en los ultimos 10 s (estabilidad)', float(t_u1[-1] - t_u1[-2]), 'K'))
        print('  Curva T(U1) vs t:', [f'{tt:.0f}s:{v:.1f}' for tt, v in zip(tiempos[::3], t_u1[::3])], flush=True)
        with (FEM / 'resultados' / 'caso_12V_2V_05A_aire.csv').open('w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['cantidad', 'valor', 'unidad'])
            for n, v, u in filas:
                w.writerow([n, f'{v:.4g}', u])
                print(f'  {n}: {v:.4g} {u}', flush=True)
        model.save(str(MODELO))
        print('Guardado:', MODELO)
        return 0
    finally:
        client.disconnect()


if __name__ == '__main__':
    raise SystemExit(main())
