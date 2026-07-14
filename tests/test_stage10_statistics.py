from syndesis.stage10.statistics import benjamini_hochberg, bootstrap_ci, paired_permutation_pvalue


def test_bootstrap_ci_and_correction():
    low, high = bootstrap_ci([0.0, 1.0, 1.0], iterations=100, seed=1)
    assert low <= high
    corrected = benjamini_hochberg([0.01, 0.04, 0.03])
    assert len(corrected) == 3
    assert all(0 <= p <= 1 for p in corrected)


def test_paired_permutation_returns_probability():
    p = paired_permutation_pvalue([1, 0, 1], [0, 0, 0])
    assert 0 <= p <= 1
