#!/usr/bin/env bash
# Everything after the workflow runs: grouping variants, scoring, figures.
# Expects the workflow outputs in work/ (run_benchmark_container.sh) and the
# truth sets in truth/ (build_truth.py).
set -euo pipefail
cd "$(dirname "$0")"
PY=${PY:-.env/bin/python}
SAMPLES=${SAMPLES:-S_coelicolor_A3_2,S_albidoflavus_J1074,S_avermitilis_MA4680,S_tropica_CNB440}
# The grouping also reads the Bakta genome FASTA (contig lengths), which work/
# does not carry; point RUNS at a work directory that has data/bakta/<sample>/<sample>.fna.
RUNS=${RUNS:-work}

for sample in ${SAMPLES//,/ }; do
  $PY run_consensus_variant.py --workdir "$RUNS" --sample "$sample" --mode voting --min-callers 2 --out "variants/$sample/voting2"
  $PY run_consensus_variant.py --workdir "$RUNS" --sample "$sample" --mode voting --min-callers 1 --out "variants/$sample/voting1"
  $PY run_consensus_variant.py --workdir "$RUNS" --sample "$sample" --mode voting --min-callers 3 --out "variants/$sample/voting3"
  $PY run_consensus_variant.py --workdir "$RUNS" --sample "$sample" --mode containment --out "variants/$sample/containment"
  # The default variant must reproduce what the workflow itself wrote.
  cmp -s "variants/$sample/voting2/consensus_bgcs.tsv" "work/results/$sample/summary/consensus_bgcs.tsv" \
    || { echo "voting2 differs from the workflow output for $sample" >&2; exit 1; }
done

$PY -m pytest -q test_validate_consensus.py
$PY validate_consensus.py --workdir work --variants variants --truth-dir truth --samples "$SAMPLES" --out metrics
$PY make_figures.py --metrics metrics --out figures
