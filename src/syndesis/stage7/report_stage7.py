from __future__ import annotations

from pathlib import Path

import pandas as pd


def write_stage7_report(paths: dict) -> Path:
    report = paths["reports"] / "07_candidate_library.html"
    master = pd.read_parquet(paths["processed"] / "candidate_library_master.parquet")
    subsets = pd.read_parquet(paths["processed"] / "screening_subsets.parquet")
    stage8 = pd.read_parquet(paths["processed"] / "stage8_screening_input.parquet")
    source_counts = master["source"].value_counts().rename_axis("source").reset_index(name="count")
    novelty = master["novelty_bucket"].value_counts().rename_axis("novelty_bucket").reset_index(name="count")
    subset_counts = subsets["screening_subset"].value_counts().rename_axis("screening_subset").reset_index(name="count")
    html = f"""<!doctype html><html><head><meta charset='utf-8'><title>Stage 7 Candidate Library</title>
<style>body{{font-family:Arial,sans-serif;margin:32px}}td,th{{border:1px solid #ccc;padding:4px 8px}}table{{border-collapse:collapse}}</style></head><body>
<h1>Stage 7 Source-aware EGFR Candidate Library</h1>
<p>Master molecules: {len(master)}. Stage 8 prepared rows: {len(stage8)}.</p>
<h2>Source Counts</h2>{source_counts.to_html(index=False)}
<h2>Novelty Buckets</h2>{novelty.to_html(index=False)}
<h2>Screening Subsets</h2>{subset_counts.to_html(index=False)}
<h2>Known Controls</h2>{master[master['screening_role'].isin(['native_pose_reference','known_activity_reference'])].head(25).to_html(index=False)}
<h2>Activity Labels</h2><p>Known activity labels are assigned only to curated activity-reference records. Vendor, generated, and manually supplied candidates remain unlabeled unless curated activity evidence is present.</p>
</body></html>"""
    report.write_text(html, encoding="utf-8")
    return report
