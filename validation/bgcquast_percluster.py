#!/usr/bin/env python3
"""Run BGC-QUAST and print the status it assigns to every reference BGC.

BGC-QUAST reports counts only. The status of each reference BGC is computed in
compare_to_ref_analyzer.compute_reference_coverage; this wrapper records what
that function returns and otherwise runs BGC-QUAST unchanged, so the per-locus
statuses are BGC-QUAST's own.

Usage: bgcquast_percluster.py OUT.json <bgc-quast arguments...>
"""
import json
import sys

from bgc_quast import compare_to_ref_analyzer
from bgc_quast.logger import Logger
from bgc_quast.main import run

STATUS = {"FULLY_RECOVERED": "full", "PARTIALLY_RECOVERED": "partial", "MISSED": "missed"}


def main():
    out = sys.argv[1]
    sys.argv = ["bgc-quast"] + sys.argv[2:]
    captured = {}
    original = compare_to_ref_analyzer.compute_reference_coverage

    def recording(*args, **kwargs):
        ref_bgcs = original(*args, **kwargs)
        for ref in ref_bgcs:
            # The reference id is "<sequence_id>_<cluster_id>"; the cluster id is the MIBiG accession.
            captured[ref.bgc_id.rsplit("_", 1)[-1]] = STATUS[ref.status.name]
        return ref_bgcs

    compare_to_ref_analyzer.compute_reference_coverage = recording
    run(Logger())
    json.dump(captured, open(out, "w"), indent=1)


if __name__ == "__main__":
    main()
