# Stage 12 Non-Claims

Mandatory wording:

- These candidates are computationally prioritized molecules, not experimentally confirmed EGFR inhibitors.
- Scores and poses are hypotheses requiring experimental validation.
- No claim is made about cellular potency, selectivity, pharmacokinetics, toxicity, or clinical relevance.
- Covalent EGFR inhibition is outside v1 scope.

Forbidden wording is enforced by `configs/stage12_candidate_dossiers.yaml` and `src/syndesis/stage12/nonclaim_generator.py`. Reports must not use terms that imply experimental inhibition, unlabeled activity, wet-lab hit validation, clinical interpretation, or experimentally unsupported selectivity.
