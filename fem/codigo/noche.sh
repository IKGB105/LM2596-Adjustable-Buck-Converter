#!/bin/bash
# Cadena de trabajos de la noche. Cada paso sigue aunque el anterior falle, y todo
# queda registrado en resultados/noche.log. Correr como servicio de systemd:
#   systemd-run --user --unit=comsol-noche -p MemoryMax=14G -p MemorySwapMax=0 bash codigo/noche.sh
cd "$(dirname "$0")/.."
PY=/home/tostada105/Proyects/Escuela/Comsol/.venv/bin/python
LOG=resultados/noche.log
paso() { echo "=== $(date '+%H:%M:%S') $1" >> "$LOG"; shift; "$@" >> "$LOG" 2>&1; echo "=== salida: $?" >> "$LOG"; }
echo "### Inicio $(date)" >> "$LOG"
paso "1. Malla del cobre solido (refinamiento automatico)" timeout 7200 "$PY" -u codigo/malla_automatica.py
paso "2. Caso 12V->2V 0.5A con cobre solido" env SUFIJO=_cobre timeout 3600 "$PY" -u codigo/caso_12V_2V_05A.py
paso "3. Construir modelo con aire simulado" env MODO=aire timeout 5400 "$PY" -u codigo/modelo_termico_step.py
paso "4. Resolver modelo con aire simulado" timeout 14400 "$PY" -u codigo/resolver_aire.py
echo "### Fin $(date)" >> "$LOG"
