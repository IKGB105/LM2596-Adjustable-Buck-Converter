"""Caso especifico: Vin = 12 V, Vout = 2 V, Iout = 0.5 A.

Carga BUCK_LM2596_termico_STEP.mph (no lo modifica), cambia solo el punto de
operacion, resuelve el estudio estacionario y guarda:
  - BUCK_LM2596_termico_12V_2V_05A.mph   (modelo resuelto de este caso)
  - resultados/caso_12V_2V_05A.csv
  - figuras/temperatura_3D_12V_2V_05A.png
Las perdidas usan el mismo modelo analitico del resto del estudio (parametros
Vsat, Vf, t_sw, DCR1, DCR2 sin cambios).
"""
from __future__ import annotations

import csv
import os
from pathlib import Path

import mph
import numpy as np

FEM = Path(__file__).resolve().parent.parent
# SUFIJO='_cobre' usa el modelo con pistas solidas (BUCK_LM2596_termico_STEP_cobre.mph)
SUFIJO = os.environ.get('SUFIJO', '')
ORIGEN = FEM / f'BUCK_LM2596_termico_STEP{SUFIJO}.mph'
DESTINO = FEM / f'BUCK_LM2596_termico_12V_2V_05A{SUFIJO}.mph'
CONDICION = {'Vin': '12[V]', 'Vout': '2[V]', 'Iout': '0.5[A]', 'con_vias': '1', 'eps_rad': '0', 'incluir_L': '1'}
TEMPS = [('T union U1 (cara tab-cuerpo)', 'max_fU1(T)'), ('T max U1', 'max_U1(T)'),
         ('T max D1', 'max_D1(T)'), ('T max placa', 'max_PCB(T)'), ('T max en todo', 'max_todo(T)')]


def main() -> int:
    client = mph.start(cores=4)
    try:
        model = client.load(str(ORIGEN))
        j = model.java
        for k, v in CONDICION.items():
            j.param().set(k, v)
        print('Parametros:', {k: str(j.param().get(k)) for k in CONDICION})
        model.solve('1 Caso nominal (3 A, con vias)')
        j.study('std1').label('Caso 12 V -> 2 V, 0.5 A (con vias, solo conveccion)')
        d = next(x for x in model / 'datasets' if x.tag() == 'dset1')
        ev = lambda e, u=None: float(np.real(np.ravel(model.evaluate(e, u, dataset=d))[0]))
        filas = [('Vin', ev('Vin', 'V'), 'V'), ('Vout', ev('Vout', 'V'), 'V'), ('Iout', ev('Iout', 'A'), 'A'),
                 ('Ciclo de trabajo D', ev('Dcic'), '')]
        filas += [(p, ev(p, 'W'), 'W') for p in ('P_U1', 'P_D1', 'P_L1', 'P_L2', 'P_tot')]
        filas += [('Potencia de salida', ev('Vout*Iout', 'W'), 'W'), ('Eficiencia estimada', 100 * ev('eficiencia'), '%')]
        filas += [(n, ev(e, 'degC'), 'degC') for n, e in TEMPS]
        (FEM / 'resultados').mkdir(exist_ok=True)
        with (FEM / 'resultados' / f'caso_12V_2V_05A{SUFIJO}.csv').open('w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['cantidad', 'valor', 'unidad'])
            for n, v, u in filas:
                w.writerow([n, f'{v:.4g}', u])
                print(f'  {n}: {v:.4g} {u}')

        if 'pg_T' in [str(t) for t in j.result().tags()]:
            j.result().remove('pg_T')
        if 'pg_T_sinvias' in [str(t) for t in j.result().tags()]:
            j.result().remove('pg_T_sinvias')
        pg = j.result().create('pg_T', 'PlotGroup3D')
        pg.label('Temperatura: 12 V -> 2 V, 0.5 A')
        pg.set('data', 'dset1')
        pg.set('titletype', 'manual')
        pg.set('title', 'Temperatura [degC]: Vin 12 V, Vout 2 V, Iout 0.5 A, con vias, solo conveccion, T amb 25 degC')
        s = pg.create('surf1', 'Surface')
        s.set('expr', 'T')
        s.set('unit', 'degC')
        s.set('colortable', 'HeatCameraLight')
        e = j.result().export().create('img_tmp', 'Image')
        e.set('plotgroup', 'pg_T')
        e.set('filename', str(FEM / 'figuras' / f'temperatura_3D_12V_2V_05A{SUFIJO}.png'))
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
        model.save(str(DESTINO))
        print('Guardado:', DESTINO)
        return 0
    finally:
        client.disconnect()


if __name__ == '__main__':
    raise SystemExit(main())
