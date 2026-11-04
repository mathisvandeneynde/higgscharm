#!/bin/bash
# use_redirector.sh — let the partitions of one or more datasets read via the CMS global redirector
# Usage:   ./use_redirector.sh <workflow> <year> <dataset> [<dataset> ...]
# Example: ./use_redirector.sh xtoemu 2022postEE TTto2L2Nu TTto4Q
#
# Every file URL in condor/<workflow>/<year>/<dataset>/partitions.json that points to a fixed
# grid site is rewritten to root://cms-xrd-global.cern.ch//store/..., so xrootd picks a site
# that has the file and responds. Files on CERN EOS (eoscms.cern.ch) are left untouched.
# A backup partitions.json.bak is kept. Afterwards resubmit with jobs_status.py:
# answer 'n' to "Regenerate" (that would rebuild partitions.json) and 'y' to "Resubmit".

set -e
WORKFLOW=$1; YEAR=$2; shift 2 || true
if [ -z "$WORKFLOW" ] || [ -z "$YEAR" ] || [ $# -eq 0 ]; then
  echo "Usage: $0 <workflow> <year> <dataset> [<dataset> ...]"; exit 1
fi

for DATASET in "$@"; do
  P="condor/${WORKFLOW}/${YEAR}/${DATASET}/partitions.json"
  if [ ! -f "$P" ]; then echo "!! $P does not exist"; continue; fi
  cp "$P" "${P}.bak"
  python3 - "$P" "$DATASET" <<'PY'
import json, re, sys
path, dataset = sys.argv[1], sys.argv[2]
REDIRECTOR = "root://cms-xrd-global.cern.ch/"
counts = {"redirector": 0, "eos": 0, "other": 0}

def fix(url):
    if "eoscms.cern.ch" in url:
        counts["eos"] += 1
        return url
    m = re.search(r"(/store/.*)$", url)
    if url.startswith("root://") and m:
        counts["redirector"] += 1
        return REDIRECTOR + m.group(1)
    counts["other"] += 1
    return url

def walk(x):
    if isinstance(x, dict):
        return {k: walk(v) for k, v in x.items()}
    if isinstance(x, list):
        return [walk(v) for v in x]
    if isinstance(x, str):
        return fix(x)
    return x

data = walk(json.load(open(path)))
json.dump(data, open(path, "w"), indent=4)
print(f"{dataset}: {counts['redirector']} via redirector, {counts['eos']} via CERN EOS, {counts['other']} unchanged")
PY
done
