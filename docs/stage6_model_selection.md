# Stage 6 Model Selection

Model selection prioritizes leakage-clean models that beat original docking score and GNINA CNNscore baselines on the primary validation split. The selected artifact stores a ranker, calibrated confidence classifier, feature schema, metrics, and model card.

If no trained model beats both baselines, Stage 6 records a failed selection reason instead of silently promoting a weaker model.
