"""Resuelve BUCK_LM2596_termico_STEP.mph (construido con modelo_termico_step.py),
exporta figuras a figuras/ y tablas a resultados/.

Correr con limite de memoria, por ejemplo:
  systemd-run --user --scope -p MemoryMax=8G python codigo/resolver_termico_step.py
"""
from __future__ import annotations

import csv
from pathlib import Path

import mph
import numpy as np

FEM = Path(__file__).resolve().parent.parent
MODELO = FEM / 'BUCK_LM2596_termico_STEP.mph'
FIG = FEM / 'figuras'
RES = FEM / 'resultados'

TEMPS = [('T union U1 (cara tab-cuerpo)', 'max_fU1(T)'), ('T max U1', 'max_U1(T)'),
         ('T max D1', 'max_D1(T)'), ('T max placa', 'max_PCB(T)'), ('T max en todo', 'max_todo(T)')]
POTENCIAS = ['P_U1', 'P_D1', 'P_L1', 'P_L2', 'P_tot']
VARIANTES = ['Con vias (base)', 'Sin vias', 'Con vias + radiacion (eps 0.9)',
             'Sin vias + radiacion (eps 0.9)', 'Con vias, sin perdidas en L1/L2']


def valor(model, expr, unidad=None, dataset=None, outer=None):
    kw = {}
    if dataset is not None:
        kw['dataset'] = dataset
    if outer is not None:
        # Un barrido de parametros que no cambian la geometria queda como soluciones
        # "internas" de un solo dataset; uno que si la cambia, como "externas"
        try:
            return float(np.real(np.ravel(model.evaluate(expr, unidad, **kw, inner=[outer]))[0]))
        except Exception:
            kw['outer'] = outer
    return float(np.real(np.ravel(model.evaluate(expr, unidad, **kw))[0]))


def grafica_T(j, tag, titulo, dset):
    if tag in [str(t) for t in j.result().tags()]:
        j.result().remove(tag)  # permite volver a correr el script
    pg = j.result().create(tag, 'PlotGroup3D')
    pg.label(titulo)
    pg.set('data', dset)
    pg.set('titletype', 'manual')
    pg.set('title', titulo)
    s = pg.create('surf1', 'Surface')
    s.set('expr', 'T')
    s.set('unit', 'degC')
    s.set('colortable', 'HeatCameraLight')
    return pg


def exportar(j, tag, archivo, vista=None):
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


def nuevo_dataset(model, antes):
    """Dataset que aparecio al resolver un barrido (el de soluciones parametricas)."""
    nuevos = [d for d in model / 'datasets' if d.tag() not in antes]
    print('  datasets nuevos:', [d.name() for d in nuevos])
    param = [d for d in nuevos if 'aram' in d.name()]
    return (param or nuevos)[-1]


def main() -> int:
    FIG.mkdir(exist_ok=True)
    RES.mkdir(exist_ok=True)
    client = mph.start(cores=4)
    try:
        model = client.load(str(MODELO))
        j = model.java

        # 1. Caso nominal
        print('Resolviendo caso nominal...', flush=True)
        model.solve('1 Caso nominal (3 A, con vias)')
        d1 = next(d for d in model / 'datasets' if 'Caso nominal' in d.name())
        pot = {p: valor(model, p, 'W', d1) for p in POTENCIAS}
        efi = valor(model, 'eficiencia', None, d1)
        tn = {n: valor(model, e, 'degC', d1) for n, e in TEMPS}
        print('Potencias [W]:', {k: round(v, 3) for k, v in pot.items()}, 'eficiencia %.1f %%' % (100 * efi))
        for n, v in tn.items():
            print(f'  {n}: {v:.1f} degC')
        with (RES / 'caso_nominal.csv').open('w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['cantidad', 'valor', 'unidad'])
            for k, v in pot.items():
                w.writerow([k, f'{v:.4f}', 'W'])
            w.writerow(['eficiencia', f'{100 * efi:.1f}', '%'])
            for n, v in tn.items():
                w.writerow([n, f'{v:.2f}', 'degC'])
        dtag = d1.tag()
        grafica_T(j, 'pg_T', 'Temperatura [degC]: Vin 35 V, Vout 1.23 V, Iout 3 A, con vias, solo conveccion', dtag)
        exportar(j, 'pg_T', FIG / 'temperatura_3D_nominal.png')
        model.save(str(MODELO))

        # 2. Barrido de corriente
        print('Resolviendo barrido de corriente...', flush=True)
        antes = {d.tag() for d in model / 'datasets'}
        model.solve('2 Barrido de corriente')
        d2 = nuevo_dataset(model, antes)
        corrientes = [0.5, 1, 1.5, 2, 2.5, 3]
        with (RES / 'barrido_corriente.csv').open('w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['Iout [A]'] + [f'{p} [W]' for p in POTENCIAS] + [f'{n} [degC]' for n, _ in TEMPS])
            for i, iout in enumerate(corrientes, 1):
                fila = [iout] + [f'{valor(model, p, "W", d2, i):.3f}' for p in POTENCIAS] \
                       + [f'{valor(model, e, "degC", d2, i):.1f}' for _, e in TEMPS]
                w.writerow(fila)
                print('  Iout', iout, 'A ->', fila[-5:])
        model.save(str(MODELO))

        # 3. Variantes
        print('Resolviendo variantes...', flush=True)
        antes = {d.tag() for d in model / 'datasets'}
        model.solve('3 Variantes: vias, radiacion, sin L')
        d3 = nuevo_dataset(model, antes)
        with (RES / 'variantes.csv').open('w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['variante'] + [f'{n} [degC]' for n, _ in TEMPS])
            for i, nombre in enumerate(VARIANTES, 1):
                fila = [nombre] + [f'{valor(model, e, "degC", d3, i):.1f}' for _, e in TEMPS]
                w.writerow(fila)
                print(' ', fila)
        pg = grafica_T(j, 'pg_T_sinvias', 'Temperatura [degC]: misma placa SIN vias (3 A)', d3.tag())
        pg.setIndex('looplevel', '2', 0)
        exportar(j, 'pg_T_sinvias', FIG / 'temperatura_3D_sin_vias.png')
        model.save(str(MODELO))
        print('Guardado:', MODELO)
        return 0
    finally:
        client.disconnect()


if __name__ == '__main__':
    raise SystemExit(main())
