"""Retrospective enrichment benchmark (LIT-PCBA EGFR) across the receptor ensemble.

This subpackage adds an external, recognized-negative validation of the Syndesis
funnel: actives vs decoys are docked across all ensemble receptors with GPU-batched
Uni-Dock, GNINA-rescored, ProLIF interaction-fingerprinted, and scored under three
matched arms (docking-score-only, +interaction-constraint, +pose-confidence). It
reuses the existing scoring/parsing code but replaces the serial, per-ligand
subprocess docking with Uni-Dock ``--gpu_batch`` (one GPU call per receptor) and
batched GNINA scoring (one container per receptor), turning a weeks-long serial CPU
run into a GPU-saturated campaign of hours.
"""
