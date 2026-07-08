# Stage 5 Interaction-Fingerprint Atlas Protocol

Stage 5 builds an EGFR interaction atlas from native co-crystal structures, computes docked-pose interaction fingerprints, measures interaction recovery, clusters binding modes, and writes final interaction-aware pose labels.

The required engine is ProLIF. If ProLIF is unavailable, Stage 5 fails instead of producing fallback interaction fingerprints.

Stage 5 does not train a pose-confidence model, run MD, generate analogs, or screen vendor candidates.
