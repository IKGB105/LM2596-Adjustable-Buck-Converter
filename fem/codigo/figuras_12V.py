"""Figuras limpias del caso 12 V -> 2 V, 0.5 A (BUCK_LM2596_termico_12V_2V_05A.mph).

Usa una vista nueva SIN planos de corte (el archivo de geometria original trae un
clip plane, "plane1", que cortaba la placa por U1 y D1 en las imagenes anteriores).
Solo exporta imagenes: no guarda cambios en el modelo.
Salidas en figuras/12V_*.png
"""
from __future__ import annotations

from pathlib import Path

import mph

FEM = Path(__file__).resolve().parent.parent
MODELO = FEM / 'BUCK_LM2596_termico_12V_2V_05A.mph'
FIG = FEM / 'figuras'
TITULO = 'Temperatura [degC]: Vin 12 V, Vout 2 V, Iout 0.5 A (0.44 W), con vias, solo conveccion, T amb 25 degC'

# nombre, camara (posicion, objetivo, arriba) en mm o None = isometrica automatica, ortografica, ajustar zoom
TOMAS = [
    ('12V_isometrica', None, False, True),
    ('12V_superior', ((50, -35, 300), (50, -35, 0), (0, 1, 0)), True, True),
    ('12V_inferior', ((50, -35, -300), (50, -35, 0), (0, 1, 0)), True, True),
    ('12V_acercamiento_U1_D1', ((5, -78, 48), (34, -34, 2), (0, 0, 1)), False, False),
]


def main() -> int:
    FIG.mkdir(exist_ok=True)
    client = mph.start(cores=2)
    try:
        model = client.load(str(MODELO))
        j = model.java
        comp = j.component('comp1')
        for tag, camara, orto, ajustar in TOMAS:
            vt = f'v_{tag}'
            v = comp.view().create(vt, 'geom1')
            for lista in ('clip', 'hideEntities', 'hideObjects'):
                for t in [str(x) for x in getattr(v, lista)().tags()]:
                    getattr(v, lista)().remove(t)
            cam = v.camera()
            if orto:
                cam.set('projection', 'orthographic')
            if camara:
                pos, obj, arriba = camara
                cam.set('position', [str(c) for c in pos])
                cam.set('target', [str(c) for c in obj])
                cam.set('up', [str(c) for c in arriba])

            pg = j.result().create(f'pg_{tag}', 'PlotGroup3D')
            pg.set('data', 'dset1')
            pg.set('view', vt)
            pg.set('titletype', 'manual')
            pg.set('title', TITULO)
            s = pg.create('surf1', 'Surface')
            s.set('expr', 'T')
            s.set('unit', 'degC')
            s.set('colortable', 'HeatCameraLight')

            e = j.result().export().create('img_tmp', 'Image')
            e.set('plotgroup', f'pg_{tag}')
            e.set('filename', str(FIG / f'{tag}.png'))
            e.set('size', 'manualweb')
            e.set('width', '1600')
            e.set('height', '1000')
            e.set('zoomextents', 'on' if ajustar else 'off')
            e.set('options3d', 'on')
            for p in ('title3d', 'legend3d', 'axisorientation'):
                e.set(p, 'on')
            e.set('grid', 'off')
            e.set('background', 'color')
            try:
                e.run()
                print('Exportada:', FIG / f'{tag}.png', flush=True)
            finally:
                j.result().export().remove('img_tmp')
        return 0
    finally:
        client.disconnect()


if __name__ == '__main__':
    raise SystemExit(main())
