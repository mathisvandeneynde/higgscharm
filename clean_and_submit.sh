#!/bin/bash
# clean_and_submit.sh — clean-slate routine for a new run
# Usage:   ./clean_and_submit.sh <workflow> <year> [<year> ...]
# Example: ./clean_and_submit.sh xtoemu 2022preEE 2022postEE
#
# Resources (override via environment variables):
#   MEMORY=16000 JOBFLAVOR=tomorrow ./clean_and_submit.sh xtoemu 2022preEE
#   JOBFLAVOR: espresso (20 min), microcentury (1 h), longlunch (2 h), workday (8 h),
#              tomorrow (1 day), testmatch (3 days), nextweek (1 week)

set -e  # stop immediately if something fails, instead of continuing with a half cleanup

WORKFLOW=$1
shift || true
YEARS=("$@")
MEMORY=${MEMORY:-12000}          # MB per job
JOBFLAVOR=${JOBFLAVOR:-workday}  # maximum wall time of a job
BASEDIR="/afs/cern.ch/user/m/mvandene/higgscharm"
TIMESTAMP=$(date +%Y%m%d_%H%M)

if [ -z "$WORKFLOW" ] || [ ${#YEARS[@]} -eq 0 ]; then
  echo "Usage: $0 <workflow> <year> [<year> ...]"
  exit 1
fi

echo "Workflow: $WORKFLOW | years: ${YEARS[*]} | memory: ${MEMORY} MB | jobflavor: ${JOBFLAVOR}"

for YEAR in "${YEARS[@]}"; do
  EOSDIR="/eos/user/m/mvandene/higgscharm/outputs/${WORKFLOW}/${YEAR}"
  LOGDIR="${BASEDIR}/condor/logs/${WORKFLOW}/${YEAR}"

  echo "=== [$YEAR] Step 1: remove only the jobs of ${WORKFLOW}/${YEAR} from the queue ==="
  # selects jobs via their log path (condor/logs/<workflow>/<year>/...), so other runs stay untouched
  condor_rm -constraint "regexp(\"/condor/logs/${WORKFLOW}/${YEAR}/\", UserLog)" 2>/dev/null || true
  sleep 2

  echo "=== [$YEAR] Step 2: rename the output directory ==="
  if [ -d "$EOSDIR" ]; then
    mv "$EOSDIR" "${EOSDIR}_OLD_${TIMESTAMP}"
    echo "Renamed to ${EOSDIR}_OLD_${TIMESTAMP}"
  else
    echo "No existing output directory found, nothing to rename."
  fi

  echo "=== [$YEAR] Step 3: clean up old condor logs ==="
  if [ -d "$LOGDIR" ]; then
    rm -rf "${LOGDIR:?}"/*
    echo "Logs cleaned in $LOGDIR"
  fi

  echo "=== [$YEAR] Step 4: fresh submission ==="
  cd "$BASEDIR"
  python3 runner.py --workflow "$WORKFLOW" --year "$YEAR" --eos \
    --memory "$MEMORY" --jobflavor "$JOBFLAVOR" --submit
done

echo "=== Queue ==="
condor_q mvandene -batch | tail -n 20