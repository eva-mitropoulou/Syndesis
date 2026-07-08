# Stage 12 Ranking And Decision Policy

Stage 12 ranks computational candidates by existing final candidate scores and preserves source labels, novelty buckets, medchem risk, pose confidence, GNINA evidence, interaction recovery, and MD status.

Allowed decision labels are defined in `src/egfr_dockingforge/stage12/schemas.py`.

Known controls remain diagnostic rows. Generated analogs rejected upstream remain negative-control rows. MD-stable labels are assigned only when Stage 11 production trajectory evidence exists.
