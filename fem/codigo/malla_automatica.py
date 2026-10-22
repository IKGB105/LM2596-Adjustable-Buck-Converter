"""Malla el modelo con cobre solido refinando solo donde COMSOL falla.

El STEP de KiCad trae detalles diminutos (esquinas donde una pista toca un pad,
anillos de via) en los que el mallador reporta "Intersecting edge elements".
Este script intenta mallar; si falla, toma las coordenadas del error, agrega un
tamano de elemento fino SOLO en una caja de +-0.6 mm alrededor y vuelve a
intentar. Asi no hace falta mallar fino toda la placa (eso no cabe en memoria).

Uso: python codigo/malla_automatica.py [modelo.mph]
"""
from __future__ import annotations

import re
import sys
import time
from pathlib import Path

import mph

FEM = Path(__file__).resolve().parent.parent
MODELO = Path(sys.argv[1]) if len(sys.argv) > 1 else FEM / 'BUCK_LM2596_termico_STEP_cobre.mph'
ESCALA_Z = 10      # mismo estiramiento en z que modelo_termico_step.py
RADIO = 0.6        # mm alrededor del error
MAX_INTENTOS = 15


def main() -> int:
    client = mph.start(cores=4)
    try:
        model = client.load(str(MODELO))
        comp = model.java.component('comp1')
        if 'mesh1' not in [str(t) for t in comp.mesh().tags()]:
            # Misma configuracion que modelo_termico_step.py
            mesh = comp.mesh().create('mesh1')
            mesh.feature('size').set('hauto', '5')
            ft = mesh.create('ftet1', 'FreeTet')
            ft.label('Tetraedros (geometria escalada en z)')
            ft.set('zscale', str(ESCALA_Z))
        mesh = comp.mesh('mesh1')
        ft = mesh.feature('ftet1')
        for intento in range(1, MAX_INTENTOS + 1):
            t0 = time.time()
            try:
                mesh.run()
                print(f'Intento {intento}: malla OK, {mesh.getNumElem()} elementos, '
                      f'calidad min {mesh.stat().getMinQuality():.3f} ({time.time() - t0:.0f} s)', flush=True)
                model.save(str(MODELO))
                print('Guardado:', MODELO)
                return 0
            except Exception as e:
                texto = str(e)
                m = re.search(r'Coordinates: \(([-\d.eE]+),\s*([-\d.eE]+),\s*([-\d.eE]+)\)', texto)
                if not m:
                    print('Error sin coordenadas, no se puede refinar:\n', texto[:800])
                    return 1
                x, y, zs = (float(v) for v in m.groups())
                z = zs / ESCALA_Z  # el error viene en la geometria escalada
                print(f'Intento {intento}: falla en ({x:.3f}, {y:.3f}, {z:.3f}) mm -> refino ahi '
                      f'({time.time() - t0:.0f} s)', flush=True)
                previos = [str(t) for t in ft.feature().tags() if str(t).startswith('ref')]
                tag = f'ref{len(previos) + 1}'  # continua la numeracion si se retoma
                b = comp.selection().create(f'sel_{tag}', 'Box')
                b.label(f'Refinamiento local {tag[3:]}')
                b.set('entitydim', '2')
                b.set('condition', 'intersects')
                for k, v in zip(('xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax'),
                                (x - RADIO, x + RADIO, y - RADIO, y + RADIO, z - 0.2, z + 0.2)):
                    b.set(k, f'{v:g}')
                sz = ft.create(tag, 'Size')
                sz.label(f'Fino alrededor de ({x:.2f}, {y:.2f})')
                sz.selection().geom('geom1', 2)
                sz.selection().named(f'sel_{tag}')
                sz.set('custom', True)
                sz.set('hmaxactive', True)
                sz.set('hmax', '0.08')
                sz.set('hminactive', True)
                sz.set('hmin', '0.004')
                model.save(str(MODELO))  # guarda el avance: si algo lo corta, retoma desde aqui
        print('No se logro mallar en', MAX_INTENTOS, 'intentos')
        model.save(str(MODELO))
        return 1
    finally:
        client.disconnect()


if __name__ == '__main__':
    raise SystemExit(main())
