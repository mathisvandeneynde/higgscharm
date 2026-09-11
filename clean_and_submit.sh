#!/bin/bash
# clean_and_submit.sh — schone-lei routine voor een nieuwe run
# Gebruik: ./clean_and_submit.sh <workflow> <year>
# Voorbeeld: ./clean_and_submit.sh ztoemu_CR 2022preEE

set -e  # stop meteen als er iets faalt, i.p.v. door te gaan met een halve cleanup

WORKFLOW=$1
YEAR=$2
BASEDIR="/afs/cern.ch/user/m/mvandene/higgscharm"
EOSDIR="/eos/user/m/mvandene/higgscharm/outputs/${WORKFLOW}/${YEAR}"
TIMESTAMP=$(date +%Y%m%d_%H%M)

if [ -z "$WORKFLOW" ] || [ -z "$YEAR" ]; then
  echo "Gebruik: $0 <workflow> <year>"
  exit 1
fi

echo "=== Stap 1: wachtrij leegmaken ==="
condor_rm mvandene 2>/dev/null || true
sleep 2
condor_q mvandene

echo "=== Stap 2: output-map hernoemen ==="
if [ -d "$EOSDIR" ]; then
  mv "$EOSDIR" "${EOSDIR}_OLD_${TIMESTAMP}"
  echo "Hernoemd naar ${EOSDIR}_OLD_${TIMESTAMP}"
else
  echo "Geen bestaande output-map gevonden, niets om te hernoemen."
fi

echo "=== Stap 3: oude condor-logs opruimen ==="
LOGDIR="${BASEDIR}/condor/logs/${WORKFLOW}/${YEAR}"
if [ -d "$LOGDIR" ]; then
  rm -rf "${LOGDIR:?}"/*
  echo "Logs opgeruimd in $LOGDIR"
fi

echo "=== Stap 4: verse submissie ==="
cd "$BASEDIR"
python3 runner.py --workflow "$WORKFLOW" --year "$YEAR" --eos --m 8000
# pas eventuele andere vaste vlaggen hier aan naar je gebruikelijke commando
