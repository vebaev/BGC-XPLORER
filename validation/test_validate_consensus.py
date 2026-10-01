"""The metric definitions in validate_consensus.py, on synthetic loci."""
import pandas as pd

from validate_consensus import score_arm

# Ten 1 kb genes back to back on one contig: g0 = 1-1000, g1 = 1001-2000, ...
GENES = pd.DataFrame({"contig": "c", "start": [i * 1000 + 1 for i in range(10)],
                      "end": [(i + 1) * 1000 for i in range(10)],
                      "locus_tag": ["g%d" % i for i in range(10)]})


def truth(*loci):
    return pd.DataFrame([{"mibig_id": m, "compound": m, "contig": "c", "start": s, "end": e,
                          "nested_in": n} for m, s, e, n in loci])


def preds(*intervals):
    return pd.DataFrame([{"contig": "c", "start": s, "end": e, "tool": "x"} for s, e in intervals])


def one(df, mibig_id="A"):
    return df.set_index("mibig_id").loc[mibig_id]


def test_exact_match_is_full_and_delimited():
    row = one(score_arm(preds((2001, 5000)), truth(("A", 2001, 5000, "")), GENES))
    assert row.category == "full" and row.delimited
    assert row.jaccard == 1.0 and row.gene_f1 == 1.0
    assert row.left_extension == 0 and row.right_extension == 0


def test_over_extension_is_full_but_loses_precision():
    row = one(score_arm(preds((1, 8000)), truth(("A", 2001, 5000, "")), GENES))
    assert row.category == "full" and not row.delimited
    assert row.left_extension == 2000 and row.right_extension == 3000
    assert row.gene_recall == 1.0 and row.gene_precision == 3 / 8


def test_three_thirds_are_split_and_full_multi():
    row = one(score_arm(preds((2001, 3000), (3001, 4000), (4001, 5000)),
                        truth(("A", 2001, 5000, "")), GENES))
    assert row.category == "full_multi" and row.split
    assert row.fragments == 3


def test_one_prediction_over_two_loci_is_fused():
    df = score_arm(preds((1, 10000)), truth(("A", 1001, 3000, ""), ("B", 6001, 9000, "")), GENES)
    assert one(df, "A").fused and one(df, "B").fused


def test_nested_locus_is_not_counted_as_fused():
    df = score_arm(preds((1001, 9000)), truth(("A", 1001, 9000, ""), ("B", 2001, 3000, "A")), GENES)
    assert not one(df, "A").fused and not one(df, "B").fused


def test_truncation_is_partial_with_negative_extension():
    row = one(score_arm(preds((2001, 3000)), truth(("A", 2001, 5000, "")), GENES))
    assert row.category == "partial"
    assert row.right_extension == -2000 and row.gene_precision == 1.0


def test_below_ten_percent_is_missed():
    row = one(score_arm(preds((1, 2200)), truth(("A", 2001, 5000, "")), GENES))
    assert row.category == "missed" and not row.detected and row.gene_f1 == 0.0
