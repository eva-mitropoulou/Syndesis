import pandas as pd
import pytest

from egfr_dockingforge.enrichment.native_prior import jaccard, native_union, parse_fingerprint, recall


def test_native_union_excludes_allosteric_receptor() -> None:
    native = pd.DataFrame(
        [
            {"receptor_id": "1m17_a_aq4_999", "fingerprint_sparse_json": '["hinge:793:HBDonor"]'},
            {"receptor_id": "6duk_c_jbj_1103", "fingerprint_sparse_json": '["allosteric:1:Hydrophobic"]'},
        ]
    )
    target, receptors = native_union(native, ["6duk_c_jbj_1103"])
    assert target == {"hinge:793:HBDonor"}
    assert receptors == ["1m17_a_aq4_999"]


def test_interaction_similarity_metrics() -> None:
    observed = parse_fingerprint('["a", "b"]')
    target = {"b", "c"}
    assert recall(observed, target) == pytest.approx(0.5)
    assert jaccard(observed, target) == pytest.approx(1 / 3)


def test_empty_native_prior_is_rejected() -> None:
    native = pd.DataFrame(
        [{"receptor_id": "6duk_c_jbj_1103", "fingerprint_sparse_json": '["a"]'}]
    )
    with pytest.raises(ValueError, match="no included complexes"):
        native_union(native, ["6duk_c_jbj_1103"])
