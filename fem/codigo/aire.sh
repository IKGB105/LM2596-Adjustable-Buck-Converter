#!/bin/bash
# Modelo con aire simulado (12 V -> 2 V, 0.5 A): construir y resolver.
#   systemd-run --user --unit=comsol-aire -p MemoryMax=14G -p MemorySwapMax=0 bash codigo/aire.sh
cd "$(dirname "$0")/.."
PY=/home/tostada105/Proyects/Escuela/Comsol/.venv/bin/python
LOG=resultados/aire.log
paso() { echo "=== $(date '+%H:%M:%S') $1" >> "$LOG"; shift; "$@" >> "$LOG" 2>&1; echo "=== salida: $?" >> "$LOG"; }
echo "### Inicio $(date)" >> "$LOG"
paso "Construir modelo con aire simulado" env MODO=aire timeout 5400 "$PY" -u codigo/modelo_termico_step.py
paso "Resolver modelo con aire simulado" timeout 14400 "$PY" -u codigo/resolver_aire.py
echo "### Fin $(date)" >> "$LOG"
