from __future__ import annotations

from pathlib import Path

import click

from egfr_dockingforge.data.stage1_structures import run_stage1
from egfr_dockingforge.stage1.cocrystal_benchmark import (
    acquire_stage1,
    build_stage1_benchmark,
    export_reference_complexes as run_export_reference_complexes,
    extract_native_ligands as run_extract_native_ligands,
    report_stage1 as run_report_stage1,
)
from egfr_dockingforge.stage2.receptor_ensemble import (
    build_receptor_features as run_build_receptor_features,
    cluster_receptor_features,
    export_receptor_ensemble as run_export_receptor_ensemble,
    report_stage2 as run_report_stage2,
    select_receptor_ensemble as run_select_receptor_ensemble,
)
from egfr_dockingforge.stage3.redocking_crossdocking import (
    build_docking_task_matrix as run_build_docking_task_matrix,
    compute_docking_rmsd as run_compute_docking_rmsd,
    label_stage3_poses as run_label_stage3_poses,
    prepare_docking_inputs as run_prepare_docking_inputs,
    report_stage3 as run_report_stage3,
    run_crossdocking as run_stage3_crossdocking,
    run_pose_sanity_checks as run_stage3_pose_sanity,
    run_redocking as run_stage3_redocking,
)
from egfr_dockingforge.stage4.rescoring import (
    build_rescoring_task_matrix as run_build_rescoring_task_matrix,
    diagnose_rescoring as run_diagnose_rescoring,
    parse_rescoring_outputs as run_parse_rescoring_outputs,
    report_stage4 as run_report_stage4,
    run_empirical_rescoring as run_stage4_empirical_rescoring,
    run_gnina_stage4 as run_stage4_gnina,
)
from egfr_dockingforge.stage5.interaction_atlas import (
    build_native_interaction_atlas as run_build_native_interaction_atlas,
    build_stage5_all as run_stage5_all,
    cluster_stage5_binding_modes as run_cluster_stage5_binding_modes,
    compute_pose_interactions as run_compute_pose_interactions,
    compute_recovery as run_compute_interaction_recovery,
    label_stage5_final_poses as run_label_stage5_final_poses,
    report_stage5 as run_report_stage5,
    run_stage5_plip_crosscheck,
)
from egfr_dockingforge.stage6.pose_model import (
    audit_pose_features_cli as run_audit_pose_features,
    build_pose_features_cli as run_build_pose_features,
    build_pose_labels_cli as run_build_pose_labels,
    calibrate_pose_confidence_cli as run_calibrate_pose_confidence,
    evaluate_pose_models_cli as run_evaluate_pose_models,
    explain_pose_model_cli as run_explain_pose_model,
    report_stage6_cli as run_report_stage6,
    run_stage6_all,
    select_pose_model_cli as run_select_pose_model,
    split_pose_data_cli as run_split_pose_data,
    train_pose_confidence_classifier_cli as run_train_pose_confidence_classifier,
    train_pose_rankers_cli as run_train_pose_rankers,
)
from egfr_dockingforge.stage7.candidate_library import (
    build_analog_series as run_stage7_build_analog_series,
    clean_activity_data as run_stage7_clean_activity_data,
    compute_candidate_similarity as run_stage7_compute_similarity,
    filter_candidate_library as run_stage7_filter_candidate_library,
    prepare_candidate_ligands as run_stage7_prepare_candidate_ligands,
    report_stage7 as run_report_stage7,
    run_stage7_all,
    select_screening_subsets as run_stage7_select_screening_subsets,
)
from egfr_dockingforge.stage7.library_export import (
    fetch_known_egfr_ligands as run_fetch_known_egfr_ligands,
    import_vendor_library as run_import_vendor_library,
    standardize_candidates as run_standardize_candidates,
)
from egfr_dockingforge.stage8.candidate_screening import (
    aggregate_candidate_scores_cli as run_aggregate_candidate_scores,
    apply_pose_confidence_cli as run_apply_pose_confidence,
    build_screening_task_matrix_cli as run_build_screening_task_matrix,
    compute_screening_interactions_cli as run_compute_screening_interactions,
    report_stage8_cli as run_report_stage8,
    rescore_screening_poses_cli as run_rescore_screening_poses,
    run_candidate_docking_cli as run_candidate_docking_stage8,
    run_stage8_all,
    select_ranked_candidates_cli as run_select_ranked_candidates,
)
from egfr_dockingforge.stage9.candidate_screening import (
    benchmark_analog_strategies_cli as run_benchmark_analog_strategies,
    detect_edit_sites_cli as run_detect_edit_sites,
    enumerate_rule_based_analogs_cli as run_enumerate_rule_based_analogs,
    report_stage9_cli as run_report_stage9,
    run_stage9_all,
    score_analog_acceptance_cli as run_score_analog_acceptance,
    screen_analog_batch_cli as run_screen_analog_batch,
    select_analog_seeds_cli as run_select_analog_seeds,
    validate_analog_batch_cli as run_validate_analog_batch,
)
from egfr_dockingforge.stage10.ablation_benchmark import (
    build_ablation_manifest_cli as run_build_ablation_manifest,
    compute_analog_benchmark_metrics_cli as run_compute_analog_benchmark_metrics,
    compute_score_hacking_metrics_cli as run_compute_score_hacking_metrics,
    make_ablation_plots_cli as run_make_ablation_plots,
    report_stage10_cli as run_report_stage10,
    run_ablation_statistics_cli as run_ablation_statistics_stage10,
    run_stage10_all,
)
from egfr_dockingforge.stage11.md_stress_test import (
    analyze_md_trajectories_cli as run_analyze_md_trajectories,
    build_md_systems_cli as run_build_md_systems,
    compute_md_interaction_persistence_cli as run_compute_md_interaction_persistence,
    parameterize_md_ligands_cli as run_parameterize_md_ligands,
    report_stage11_cli as run_report_stage11,
    run_md_equilibration_cli as run_stage11_md_equilibration,
    run_md_minimization_cli as run_stage11_md_minimization,
    run_md_production_cli as run_stage11_md_production,
    run_stage11_all,
    score_md_stability_cli as run_score_md_stability,
    select_md_candidates_cli as run_select_md_candidates,
)
from egfr_dockingforge.stage12.candidate_card_builder import build_candidate_cards as run_build_candidate_cards
from egfr_dockingforge.stage12.candidate_dossiers import run_stage12_all
from egfr_dockingforge.stage12.dataset_card import build_dataset_card as run_build_dataset_card
from egfr_dockingforge.stage12.dossier_renderer import render_candidate_dossiers as run_render_candidate_dossiers
from egfr_dockingforge.stage12.final_candidate_selection import build_final_candidate_table as run_build_final_candidate_table
from egfr_dockingforge.stage12.model_card import build_model_card as run_build_model_card
from egfr_dockingforge.stage12.provenance_bundle import build_provenance_bundle as run_build_provenance_bundle
from egfr_dockingforge.stage12.report_stage12 import render_final_figures as run_render_final_figures
from egfr_dockingforge.stage12.report_stage12 import report_stage12 as run_report_stage12
from egfr_dockingforge.stage0.scope_schema import validate_scope_files


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
def main() -> None:
    """EGFR DockingForge command-line interface."""


@main.command("validate-scope")
@click.option("--scope", "scope_path", default="configs/project_scope.yaml", show_default=True)
@click.option("--sources", "sources_path", default="data/references/stage0_sources.yaml", show_default=True)
def validate_scope(scope_path: str, sources_path: str) -> None:
    """Validate Stage 0 scope decisions and source traceability."""
    result = validate_scope_files(scope_path, sources_path)
    if not result.valid:
        message = "\n".join(f"- {error}" for error in result.errors)
        raise click.ClickException(f"Stage 0 scope validation failed:\n{message}")
    click.echo("Stage 0 scope validation passed.")


@main.command("fetch-egfr-cocrystals")
@click.option("--config", "config_path", default="configs/stage1_cocrystal_benchmark.yaml", show_default=True)
def fetch_egfr_cocrystals(config_path: str) -> None:
    """Fetch Stage 1 candidate EGFR co-crystal source files."""
    summary = acquire_stage1(Path(config_path))
    click.echo(f"Fetched/acquired {len(summary['pdb_ids'])} candidate PDB entries.")
    click.echo(f"Manifest: {summary['manifest_path']}")


@main.command("build-cocrystal-benchmark")
@click.option("--config", "config_path", default="configs/stage1_cocrystal_benchmark.yaml", show_default=True)
def build_cocrystal_benchmark(config_path: str) -> None:
    """Build the curated Stage 1 EGFR co-crystal benchmark."""
    summary = build_stage1_benchmark(Path(config_path))
    click.echo(f"Stage 1 status: {summary['status']}")
    click.echo(f"Candidate complexes: {summary['candidate_complexes']}")
    click.echo(f"Retained complexes: {summary['retained_complexes']}")
    click.echo(f"Rejected complexes: {summary['rejected_complexes']}")
    click.echo(f"Report: {summary['report']}")


@main.command("extract-native-ligands")
@click.option("--config", "config_path", default="configs/stage1_cocrystal_benchmark.yaml", show_default=True)
def extract_native_ligands(config_path: str) -> None:
    """Extract native ligand coordinate files for Stage 1 complexes."""
    summary = run_extract_native_ligands(Path(config_path))
    click.echo(f"Native ligand extraction complete: {summary['candidate_complexes']} complexes processed.")


@main.command("export-reference-complexes")
@click.option("--config", "config_path", default="configs/stage1_cocrystal_benchmark.yaml", show_default=True)
def export_reference_complexes(config_path: str) -> None:
    """Export native complexes, clean receptors, ligands, and pocket-water files."""
    summary = run_export_reference_complexes(Path(config_path))
    click.echo(f"Reference complex export complete: {summary['candidate_complexes']} complexes processed.")


@main.command("report-stage1")
@click.option("--config", "config_path", default="configs/stage1_cocrystal_benchmark.yaml", show_default=True)
def report_stage1(config_path: str) -> None:
    """Generate the Stage 1 cocrystal benchmark HTML report."""
    summary = run_report_stage1(Path(config_path))
    click.echo(f"Stage 1 report: {summary['report']}")


@main.command("build-receptor-features")
@click.option("--config", "config_path", default="configs/stage2_receptor_ensemble.yaml", show_default=True)
def build_receptor_features(config_path: str) -> None:
    """Build Stage 2 receptor feature tables."""
    summary = run_build_receptor_features(Path(config_path))
    click.echo(f"Stage 2 receptor features: {summary['receptor_count']} receptors.")


@main.command("cluster-receptors")
@click.option("--config", "config_path", default="configs/stage2_receptor_ensemble.yaml", show_default=True)
def cluster_receptors(config_path: str) -> None:
    """Cluster Stage 2 receptors within state strata."""
    summary = cluster_receptor_features(Path(config_path))
    click.echo(f"Stage 2 receptor clusters: {summary['cluster_count']} clusters.")


@main.command("select-receptor-ensemble")
@click.option("--config", "config_path", default="configs/stage2_receptor_ensemble.yaml", show_default=True)
def select_receptor_ensemble(config_path: str) -> None:
    """Select the Stage 2 receptor ensemble v1."""
    summary = run_select_receptor_ensemble(Path(config_path))
    click.echo(f"Selected receptors: {summary['selected_receptors']}")
    click.echo(f"Holdout receptors: {summary['holdout_receptors']}")


@main.command("export-receptor-ensemble")
@click.option("--config", "config_path", default="configs/stage2_receptor_ensemble.yaml", show_default=True)
def export_receptor_ensemble(config_path: str) -> None:
    """Export selected Stage 2 receptor ensemble files."""
    summary = run_export_receptor_ensemble(Path(config_path))
    click.echo(f"Exported selected receptors: {summary['selected_receptors']}")


@main.command("report-stage2")
@click.option("--config", "config_path", default="configs/stage2_receptor_ensemble.yaml", show_default=True)
def report_stage2(config_path: str) -> None:
    """Generate the Stage 2 receptor ensemble report."""
    summary = run_report_stage2(Path(config_path))
    click.echo(f"Stage 2 report: {summary['report']}")


@main.command("prepare-docking-inputs")
@click.option("--config", "config_path", default="configs/stage3_redocking_crossdocking.yaml", show_default=True)
def prepare_docking_inputs(config_path: str) -> None:
    summary = run_prepare_docking_inputs(Path(config_path))
    click.echo(f"Prepared receptors: {summary['receptors']}; ligands: {summary['ligands']}")


@main.command("build-docking-task-matrix")
@click.option("--config", "config_path", default="configs/stage3_redocking_crossdocking.yaml", show_default=True)
def build_docking_task_matrix(config_path: str) -> None:
    summary = run_build_docking_task_matrix(Path(config_path))
    click.echo(f"Docking tasks: {summary['tasks']}")


@main.command("run-redocking")
@click.option("--config", "config_path", default="configs/stage3_redocking_crossdocking.yaml", show_default=True)
def run_redocking(config_path: str) -> None:
    summary = run_stage3_redocking(Path(config_path))
    click.echo(f"Docking runs recorded: {summary['runs']}; poses parsed: {summary['poses']}")


@main.command("run-crossdocking")
@click.option("--config", "config_path", default="configs/stage3_redocking_crossdocking.yaml", show_default=True)
def run_crossdocking(config_path: str) -> None:
    summary = run_stage3_crossdocking(Path(config_path))
    click.echo(f"Docking runs recorded: {summary['runs']}; poses parsed: {summary['poses']}")


@main.command("compute-docking-rmsd")
@click.option("--config", "config_path", default="configs/stage3_redocking_crossdocking.yaml", show_default=True)
def compute_docking_rmsd(config_path: str) -> None:
    summary = run_compute_docking_rmsd(Path(config_path))
    click.echo(f"RMSD rows: {summary['poses']}")


@main.command("run-pose-sanity-checks")
@click.option("--config", "config_path", default="configs/stage3_redocking_crossdocking.yaml", show_default=True)
def run_pose_sanity_checks(config_path: str) -> None:
    summary = run_stage3_pose_sanity(Path(config_path))
    click.echo(f"Sanity rows: {summary['poses']}")


@main.command("label-stage3-poses")
@click.option("--config", "config_path", default="configs/stage3_redocking_crossdocking.yaml", show_default=True)
def label_stage3_poses(config_path: str) -> None:
    summary = run_label_stage3_poses(Path(config_path))
    click.echo(f"Pose labels: {summary['labels']}; task metrics: {summary['task_metrics']}")


@main.command("report-stage3")
@click.option("--config", "config_path", default="configs/stage3_redocking_crossdocking.yaml", show_default=True)
def report_stage3(config_path: str) -> None:
    summary = run_report_stage3(Path(config_path))
    click.echo(f"Stage 3 report: {summary['report']}")


@main.command("build-rescoring-task-matrix")
@click.option("--config", "config_path", default="configs/stage4_rescoring.yaml", show_default=True)
def build_rescoring_task_matrix(config_path: str) -> None:
    summary = run_build_rescoring_task_matrix(Path(config_path))
    click.echo(f"Rescoring tasks: {summary['tasks']}; ready: {summary['ready']}")


@main.command("run-gnina-rescoring")
@click.option("--config", "config_path", default="configs/stage4_rescoring.yaml", show_default=True)
def run_gnina_rescoring(config_path: str) -> None:
    summary = run_stage4_gnina(Path(config_path))
    click.echo(f"GNINA runs: {summary['runs']}; scores: {summary['scores']}")


@main.command("run-empirical-rescoring")
@click.option("--config", "config_path", default="configs/stage4_rescoring.yaml", show_default=True)
def run_empirical_rescoring(config_path: str) -> None:
    summary = run_stage4_empirical_rescoring(Path(config_path))
    click.echo(f"Empirical score rows: {summary['scores']}")


@main.command("parse-rescoring-outputs")
@click.option("--config", "config_path", default="configs/stage4_rescoring.yaml", show_default=True)
def parse_rescoring_outputs(config_path: str) -> None:
    summary = run_parse_rescoring_outputs(Path(config_path))
    click.echo(f"Pose score rows: {summary['pose_scores']}")


@main.command("diagnose-rescoring")
@click.option("--config", "config_path", default="configs/stage4_rescoring.yaml", show_default=True)
def diagnose_rescoring(config_path: str) -> None:
    summary = run_diagnose_rescoring(Path(config_path))
    click.echo(f"Rescoring task metrics: {summary['task_metrics']}; failures: {summary['failures']}")


@main.command("report-stage4")
@click.option("--config", "config_path", default="configs/stage4_rescoring.yaml", show_default=True)
def report_stage4(config_path: str) -> None:
    summary = run_report_stage4(Path(config_path))
    click.echo(f"Stage 4 report: {summary['report']}")


@main.command("build-native-interaction-atlas")
@click.option("--config", "config_path", default="configs/stage5_interaction_atlas.yaml", show_default=True)
def build_native_interaction_atlas(config_path: str) -> None:
    summary = run_build_native_interaction_atlas(Path(config_path))
    click.echo(f"Native interactions: {summary['native_interactions']}; key interactions: {summary['key_interactions']}")


@main.command("compute-pose-interactions")
@click.option("--config", "config_path", default="configs/stage5_interaction_atlas.yaml", show_default=True)
def compute_pose_interactions(config_path: str) -> None:
    summary = run_compute_pose_interactions(Path(config_path))
    click.echo(f"Pose interactions: {summary['pose_interactions']}; pose fingerprints: {summary['pose_fingerprints']}")


@main.command("compute-interaction-recovery")
@click.option("--config", "config_path", default="configs/stage5_interaction_atlas.yaml", show_default=True)
def compute_interaction_recovery(config_path: str) -> None:
    summary = run_compute_interaction_recovery(Path(config_path))
    click.echo(f"Interaction recovery rows: {summary['recovery_rows']}")


@main.command("cluster-binding-modes")
@click.option("--config", "config_path", default="configs/stage5_interaction_atlas.yaml", show_default=True)
def cluster_binding_modes(config_path: str) -> None:
    summary = run_cluster_stage5_binding_modes(Path(config_path))
    click.echo(f"Binding-mode clusters: {summary['clusters']}; rows: {summary['rows']}")


@main.command("label-final-poses")
@click.option("--config", "config_path", default="configs/stage5_interaction_atlas.yaml", show_default=True)
def label_final_poses(config_path: str) -> None:
    summary = run_label_stage5_final_poses(Path(config_path))
    click.echo(f"Final pose labels: {summary['final_labels']}; Stage 6 features: {summary['stage6_features']}")


@main.command("run-plip-crosscheck")
@click.option("--config", "config_path", default="configs/stage5_interaction_atlas.yaml", show_default=True)
def run_plip_crosscheck(config_path: str) -> None:
    summary = run_stage5_plip_crosscheck(Path(config_path))
    click.echo(f"PLIP cross-check rows: {summary['plip_rows']}")


@main.command("report-stage5")
@click.option("--config", "config_path", default="configs/stage5_interaction_atlas.yaml", show_default=True)
def report_stage5(config_path: str) -> None:
    summary = run_report_stage5(Path(config_path))
    click.echo(f"Stage 5 report: {summary['report']}")


@main.command("run-stage5")
@click.option("--config", "config_path", default="configs/stage5_interaction_atlas.yaml", show_default=True)
def run_stage5(config_path: str) -> None:
    summary = run_stage5_all(Path(config_path))
    click.echo(f"Stage 5 complete: {summary['final_labels']} final labels; report: {summary['report']}")


@main.command("build-pose-features")
@click.option("--config", "config_path", default="configs/stage6_pose_model.yaml", show_default=True)
def build_pose_features(config_path: str) -> None:
    summary = run_build_pose_features(Path(config_path))
    click.echo(f"Stage 6 feature rows: {summary['rows']}; columns: {summary['columns']}")


@main.command("audit-pose-features")
@click.option("--config", "config_path", default="configs/stage6_pose_model.yaml", show_default=True)
def audit_pose_features(config_path: str) -> None:
    summary = run_audit_pose_features(Path(config_path))
    click.echo(f"Stage 6 leakage audit features: {summary['features']}; trainable: {summary['trainable']}")


@main.command("build-pose-labels")
@click.option("--config", "config_path", default="configs/stage6_pose_model.yaml", show_default=True)
def build_pose_labels(config_path: str) -> None:
    summary = run_build_pose_labels(Path(config_path))
    click.echo(f"Stage 6 pose labels: {summary['labels']}")


@main.command("split-pose-data")
@click.option("--config", "config_path", default="configs/stage6_pose_model.yaml", show_default=True)
def split_pose_data(config_path: str) -> None:
    summary = run_split_pose_data(Path(config_path))
    click.echo(f"Stage 6 ranking groups: {summary['groups']}; split rows: {summary['split_rows']}")


@main.command("train-pose-rankers")
@click.option("--config", "config_path", default="configs/stage6_pose_model.yaml", show_default=True)
def train_pose_rankers(config_path: str) -> None:
    summary = run_train_pose_rankers(Path(config_path))
    click.echo(f"Stage 6 ranker artifact: {summary['artifact']}")


@main.command("train-pose-confidence-classifier")
@click.option("--config", "config_path", default="configs/stage6_pose_model.yaml", show_default=True)
def train_pose_confidence_classifier(config_path: str) -> None:
    summary = run_train_pose_confidence_classifier(Path(config_path))
    click.echo(f"Stage 6 confidence artifact: {summary['artifact']}")


@main.command("calibrate-pose-confidence")
@click.option("--config", "config_path", default="configs/stage6_pose_model.yaml", show_default=True)
def calibrate_pose_confidence(config_path: str) -> None:
    summary = run_calibrate_pose_confidence(Path(config_path))
    click.echo(f"Stage 6 calibration rows: {summary['calibration_rows']}")


@main.command("evaluate-pose-models")
@click.option("--config", "config_path", default="configs/stage6_pose_model.yaml", show_default=True)
def evaluate_pose_models(config_path: str) -> None:
    summary = run_evaluate_pose_models(Path(config_path))
    click.echo(f"Stage 6 metric rows: {summary['metrics']}")


@main.command("explain-pose-model")
@click.option("--config", "config_path", default="configs/stage6_pose_model.yaml", show_default=True)
def explain_pose_model(config_path: str) -> None:
    summary = run_explain_pose_model(Path(config_path))
    click.echo(f"Stage 6 feature importance rows: {summary['feature_importance_rows']}; ablation rows: {summary['ablation_rows']}")


@main.command("select-pose-model")
@click.option("--config", "config_path", default="configs/stage6_pose_model.yaml", show_default=True)
def select_pose_model(config_path: str) -> None:
    summary = run_select_pose_model(Path(config_path))
    click.echo(f"Stage 6 selected ranker: {summary['selected_ranker']}")


@main.command("report-stage6")
@click.option("--config", "config_path", default="configs/stage6_pose_model.yaml", show_default=True)
def report_stage6(config_path: str) -> None:
    summary = run_report_stage6(Path(config_path))
    click.echo(f"Stage 6 report: {summary['report']}")


@main.command("run-stage6")
@click.option("--config", "config_path", default="configs/stage6_pose_model.yaml", show_default=True)
def run_stage6(config_path: str) -> None:
    summary = run_stage6_all(Path(config_path))
    click.echo(f"Stage 6 complete: {summary['features']} feature rows; report: {summary['report']}")


@main.command("fetch-known-egfr-ligands")
@click.option("--config", "config_path", default="configs/stage7_candidate_library.yaml", show_default=True)
def fetch_known_egfr_ligands(config_path: str) -> None:
    summary = run_fetch_known_egfr_ligands(Path(config_path))
    click.echo(f"Stage 7 known ligand raw rows: ChEMBL {summary['chembl_rows']}; BindingDB {summary['bindingdb_rows']}; native {summary['native_rows']}")


@main.command("import-vendor-library")
@click.option("--config", "config_path", default="configs/stage7_candidate_library.yaml", show_default=True)
def import_vendor_library(config_path: str) -> None:
    summary = run_import_vendor_library(Path(config_path))
    click.echo(f"Stage 7 vendor rows: {summary['vendor_rows']}; manifest rows: {summary['manifest_rows']}")


@main.command("import-analogs")
@click.option("--config", "config_path", default="configs/stage7_candidate_library.yaml", show_default=True)
def import_analogs(config_path: str) -> None:
    click.echo("Stage 7 analog import complete: no generated/manual analog files configured.")


@main.command("standardize-candidates")
@click.option("--config", "config_path", default="configs/stage7_candidate_library.yaml", show_default=True)
def standardize_candidates(config_path: str) -> None:
    summary = run_standardize_candidates(Path(config_path))
    click.echo(f"Stage 7 standardized rows: {summary['standardized_rows']}; measurements: {summary['measurement_rows']}")


@main.command("clean-egfr-activity-data")
@click.option("--config", "config_path", default="configs/stage7_candidate_library.yaml", show_default=True)
def clean_egfr_activity_data(config_path: str) -> None:
    summary = run_stage7_clean_activity_data(Path(config_path))
    click.echo(f"Stage 7 activity measurements: {summary['measurement_rows']}")


@main.command("build-analog-series")
@click.option("--config", "config_path", default="configs/stage7_candidate_library.yaml", show_default=True)
def build_analog_series(config_path: str) -> None:
    summary = run_stage7_build_analog_series(Path(config_path))
    click.echo(f"Stage 7 analog series rows: {summary['series_rows']}")


@main.command("compute-candidate-similarity")
@click.option("--config", "config_path", default="configs/stage7_candidate_library.yaml", show_default=True)
def compute_candidate_similarity(config_path: str) -> None:
    summary = run_stage7_compute_similarity(Path(config_path))
    click.echo(f"Stage 7 master rows: {summary['master_rows']}; Stage 8 rows: {summary['stage8_rows']}")


@main.command("filter-candidate-library")
@click.option("--config", "config_path", default="configs/stage7_candidate_library.yaml", show_default=True)
def filter_candidate_library(config_path: str) -> None:
    summary = run_stage7_filter_candidate_library(Path(config_path))
    click.echo(f"Stage 7 filtered master rows: {summary['master_rows']}")


@main.command("select-screening-subsets")
@click.option("--config", "config_path", default="configs/stage7_candidate_library.yaml", show_default=True)
def select_screening_subsets(config_path: str) -> None:
    summary = run_stage7_select_screening_subsets(Path(config_path))
    click.echo(f"Stage 7 subset Stage 8 rows: {summary['stage8_rows']}")


@main.command("prepare-candidate-ligands")
@click.option("--config", "config_path", default="configs/stage7_candidate_library.yaml", show_default=True)
def prepare_candidate_ligands(config_path: str) -> None:
    summary = run_stage7_prepare_candidate_ligands(Path(config_path))
    click.echo(f"Stage 7 prepared Stage 8 rows: {summary['stage8_rows']}")


@main.command("report-stage7")
@click.option("--config", "config_path", default="configs/stage7_candidate_library.yaml", show_default=True)
def report_stage7(config_path: str) -> None:
    summary = run_report_stage7(Path(config_path))
    click.echo(f"Stage 7 report: {summary['report']}")


@main.command("run-stage7")
@click.option("--config", "config_path", default="configs/stage7_candidate_library.yaml", show_default=True)
def run_stage7(config_path: str) -> None:
    summary = run_stage7_all(Path(config_path))
    click.echo(f"Stage 7 complete: {summary['report']}")


@main.command("build-screening-task-matrix")
@click.option("--config", "config_path", default="configs/stage8_candidate_screening.yaml", show_default=True)
def build_screening_task_matrix(config_path: str) -> None:
    summary = run_build_screening_task_matrix(Path(config_path))
    click.echo(f"Stage 8 manifest rows: {summary['manifest_rows']}; tasks: {summary['tasks']}")


@main.command("run-candidate-docking")
@click.option("--config", "config_path", default="configs/stage8_candidate_screening.yaml", show_default=True)
def run_candidate_docking(config_path: str) -> None:
    summary = run_candidate_docking_stage8(Path(config_path))
    click.echo(f"Stage 8 docking runs: {summary['runs']}; poses: {summary['poses']}")


@main.command("rescore-screening-poses")
@click.option("--config", "config_path", default="configs/stage8_candidate_screening.yaml", show_default=True)
def rescore_screening_poses(config_path: str) -> None:
    summary = run_rescore_screening_poses(Path(config_path))
    click.echo(f"Stage 8 GNINA score rows: {summary['scores']}")


@main.command("compute-screening-interactions")
@click.option("--config", "config_path", default="configs/stage8_candidate_screening.yaml", show_default=True)
def compute_screening_interactions(config_path: str) -> None:
    summary = run_compute_screening_interactions(Path(config_path))
    click.echo(f"Stage 8 interactions: {summary['interactions']}; feature rows: {summary['features']}")


@main.command("apply-pose-confidence")
@click.option("--config", "config_path", default="configs/stage8_candidate_screening.yaml", show_default=True)
def apply_pose_confidence(config_path: str) -> None:
    summary = run_apply_pose_confidence(Path(config_path))
    click.echo(f"Stage 8 confidence rows: {summary['confidence_rows']}")


@main.command("aggregate-candidate-scores")
@click.option("--config", "config_path", default="configs/stage8_candidate_screening.yaml", show_default=True)
def aggregate_candidate_scores(config_path: str) -> None:
    summary = run_aggregate_candidate_scores(Path(config_path))
    click.echo(f"Stage 8 aggregate rows: {summary['aggregate_rows']}")


@main.command("select-ranked-candidates")
@click.option("--config", "config_path", default="configs/stage8_candidate_screening.yaml", show_default=True)
def select_ranked_candidates(config_path: str) -> None:
    summary = run_select_ranked_candidates(Path(config_path))
    click.echo(f"Stage 8 ranked rows: {summary['ranked']}; diagnostics: {summary['diagnostics']}")


@main.command("report-stage8")
@click.option("--config", "config_path", default="configs/stage8_candidate_screening.yaml", show_default=True)
def report_stage8(config_path: str) -> None:
    summary = run_report_stage8(Path(config_path))
    click.echo(f"Stage 8 report: {summary['report']}")


@main.command("run-stage8")
@click.option("--config", "config_path", default="configs/stage8_candidate_screening.yaml", show_default=True)
def run_stage8(config_path: str) -> None:
    summary = run_stage8_all(Path(config_path))
    click.echo(f"Stage 8 complete: {summary['report']}")


@main.command("select-analog-seeds")
@click.option("--config", "config_path", default="configs/stage9_deterministic_analogs.yaml", show_default=True)
def select_analog_seeds(config_path: str) -> None:
    summary = run_select_analog_seeds(Path(config_path))
    click.echo(f"Stage 9 selected seeds: {summary['seeds']}")


@main.command("detect-edit-sites")
@click.option("--config", "config_path", default="configs/stage9_deterministic_analogs.yaml", show_default=True)
def detect_edit_sites(config_path: str) -> None:
    summary = run_detect_edit_sites(Path(config_path))
    click.echo(f"Stage 9 edit sites: {summary['edit_sites']}")


@main.command("enumerate-rule-based-analogs")
@click.option("--config", "config_path", default="configs/stage9_deterministic_analogs.yaml", show_default=True)
def enumerate_rule_based_analogs(config_path: str) -> None:
    summary = run_enumerate_rule_based_analogs(Path(config_path))
    click.echo(f"Stage 9 rule-based analogs: {summary['analogs']}")


@main.command("validate-analog-batch")
@click.option("--config", "config_path", default="configs/stage9_deterministic_analogs.yaml", show_default=True)
def validate_analog_batch(config_path: str) -> None:
    summary = run_validate_analog_batch(Path(config_path))
    click.echo(f"Stage 9 validation rows: {summary['validation_rows']}")


@main.command("screen-analog-batch")
@click.option("--config", "config_path", default="configs/stage9_deterministic_analogs.yaml", show_default=True)
def screen_analog_batch(config_path: str) -> None:
    summary = run_screen_analog_batch(Path(config_path))
    click.echo(f"Stage 9 screened analogs: {summary['screened']}")


@main.command("score-analog-acceptance")
@click.option("--config", "config_path", default="configs/stage9_deterministic_analogs.yaml", show_default=True)
def score_analog_acceptance(config_path: str) -> None:
    summary = run_score_analog_acceptance(Path(config_path))
    click.echo(f"Stage 9 acceptance rows: {summary['acceptance_rows']}; accepted: {summary['accepted']}")


@main.command("benchmark-analog-strategies")
@click.option("--config", "config_path", default="configs/stage9_deterministic_analogs.yaml", show_default=True)
def benchmark_analog_strategies(config_path: str) -> None:
    summary = run_benchmark_analog_strategies(Path(config_path))
    click.echo(f"Stage 9 benchmark strategies: {summary['strategies']}")


@main.command("report-stage9")
@click.option("--config", "config_path", default="configs/stage9_deterministic_analogs.yaml", show_default=True)
def report_stage9(config_path: str) -> None:
    summary = run_report_stage9(Path(config_path))
    click.echo(f"Stage 9 report: {summary['report']}")


@main.command("run-stage9")
@click.option("--config", "config_path", default="configs/stage9_deterministic_analogs.yaml", show_default=True)
def run_stage9(config_path: str) -> None:
    summary = run_stage9_all(Path(config_path))
    click.echo(f"Stage 9 complete: {summary['report']}")


@main.command("build-ablation-manifest")
@click.option("--config", "config_path", default="configs/stage10_ablation_benchmark.yaml", show_default=True)
def build_ablation_manifest(config_path: str) -> None:
    summary = run_build_ablation_manifest(Path(config_path))
    click.echo(f"Stage 10 strategies: {summary['strategies']}; budget rows: {summary['budget_rows']}")


@main.command("compute-analog-benchmark-metrics")
@click.option("--config", "config_path", default="configs/stage10_ablation_benchmark.yaml", show_default=True)
def compute_analog_benchmark_metrics(config_path: str) -> None:
    summary = run_compute_analog_benchmark_metrics(Path(config_path))
    click.echo(f"Stage 10 metrics: seed rows {summary['seed_rows']}; strategy rows {summary['strategy_rows']}")


@main.command("compute-score-hacking-metrics")
@click.option("--config", "config_path", default="configs/stage10_ablation_benchmark.yaml", show_default=True)
def compute_score_hacking_metrics(config_path: str) -> None:
    summary = run_compute_score_hacking_metrics(Path(config_path))
    click.echo(f"Stage 10 score-hacking cases: {summary['cases']}")


@main.command("run-ablation-statistics")
@click.option("--config", "config_path", default="configs/stage10_ablation_benchmark.yaml", show_default=True)
def run_ablation_statistics(config_path: str) -> None:
    summary = run_ablation_statistics_stage10(Path(config_path))
    click.echo(f"Stage 10 comparisons: {summary['comparisons']}; ablations: {summary['ablations']}")


@main.command("make-ablation-plots")
@click.option("--config", "config_path", default="configs/stage10_ablation_benchmark.yaml", show_default=True)
def make_ablation_plots(config_path: str) -> None:
    summary = run_make_ablation_plots(Path(config_path))
    click.echo(f"Stage 10 figures: {summary['figures']}")


@main.command("report-stage10")
@click.option("--config", "config_path", default="configs/stage10_ablation_benchmark.yaml", show_default=True)
def report_stage10(config_path: str) -> None:
    summary = run_report_stage10(Path(config_path))
    click.echo(f"Stage 10 report: {summary['report']}")


@main.command("run-stage10")
@click.option("--config", "config_path", default="configs/stage10_ablation_benchmark.yaml", show_default=True)
def run_stage10(config_path: str) -> None:
    summary = run_stage10_all(Path(config_path))
    click.echo(f"Stage 10 complete: {summary['report']}")


@main.command("select-md-candidates")
@click.option("--config", "config_path", default="configs/stage11_md_stress_test.yaml", show_default=True)
def select_md_candidates(config_path: str) -> None:
    summary = run_select_md_candidates(Path(config_path))
    click.echo(f"Stage 11 MD candidates: {summary['candidates']}")


@main.command("parameterize-md-ligands")
@click.option("--config", "config_path", default="configs/stage11_md_stress_test.yaml", show_default=True)
def parameterize_md_ligands(config_path: str) -> None:
    summary = run_parameterize_md_ligands(Path(config_path))
    click.echo(f"Stage 11 ligand parameterization ready: {summary['parameterized']}; failed: {summary['failed']}")


@main.command("build-md-systems")
@click.option("--config", "config_path", default="configs/stage11_md_stress_test.yaml", show_default=True)
def build_md_systems(config_path: str) -> None:
    summary = run_build_md_systems(Path(config_path))
    click.echo(f"Stage 11 systems: {summary['systems']}; ready: {summary['ready']}")


@main.command("run-md-minimization")
@click.option("--config", "config_path", default="configs/stage11_md_stress_test.yaml", show_default=True)
def run_md_minimization(config_path: str) -> None:
    summary = run_stage11_md_minimization(Path(config_path))
    click.echo(f"Stage 11 MD run rows: {summary['run_rows']}")


@main.command("run-md-equilibration")
@click.option("--config", "config_path", default="configs/stage11_md_stress_test.yaml", show_default=True)
def run_md_equilibration(config_path: str) -> None:
    summary = run_stage11_md_equilibration(Path(config_path))
    click.echo(f"Stage 11 MD run rows: {summary['run_rows']}")


@main.command("run-md-production")
@click.option("--config", "config_path", default="configs/stage11_md_stress_test.yaml", show_default=True)
def run_md_production(config_path: str) -> None:
    summary = run_stage11_md_production(Path(config_path))
    click.echo(f"Stage 11 MD run rows: {summary['run_rows']}")


@main.command("analyze-md-trajectories")
@click.option("--config", "config_path", default="configs/stage11_md_stress_test.yaml", show_default=True)
def analyze_md_trajectories(config_path: str) -> None:
    summary = run_analyze_md_trajectories(Path(config_path))
    click.echo(f"Stage 11 MD metrics: {summary['metrics']}; QC rows: {summary['qc']}")


@main.command("compute-md-interaction-persistence")
@click.option("--config", "config_path", default="configs/stage11_md_stress_test.yaml", show_default=True)
def compute_md_interaction_persistence(config_path: str) -> None:
    summary = run_compute_md_interaction_persistence(Path(config_path))
    click.echo(f"Stage 11 MD interactions: {summary['interactions']}; summary: {summary['summary']}")


@main.command("score-md-stability")
@click.option("--config", "config_path", default="configs/stage11_md_stress_test.yaml", show_default=True)
def score_md_stability(config_path: str) -> None:
    summary = run_score_md_stability(Path(config_path))
    click.echo(f"Stage 11 stability labels: {summary['labels']}; Stage 10B rows: {summary['post_md_rows']}")


@main.command("report-stage11")
@click.option("--config", "config_path", default="configs/stage11_md_stress_test.yaml", show_default=True)
def report_stage11(config_path: str) -> None:
    summary = run_report_stage11(Path(config_path))
    click.echo(f"Stage 11 report: {summary['report']}")


@main.command("run-stage11")
@click.option("--config", "config_path", default="configs/stage11_md_stress_test.yaml", show_default=True)
def run_stage11(config_path: str) -> None:
    summary = run_stage11_all(Path(config_path))
    click.echo(f"Stage 11 complete: {summary['report']}")


@main.command("build-final-candidate-table")
@click.option("--config", "config_path", default="configs/stage12_candidate_dossiers.yaml", show_default=True)
def build_final_candidate_table(config_path: str) -> None:
    summary = run_build_final_candidate_table(Path(config_path))
    click.echo(f"Stage 12 final candidates: {summary['candidates']}; selected dossiers: {summary['selected']}")


@main.command("build-candidate-cards")
@click.option("--config", "config_path", default="configs/stage12_candidate_dossiers.yaml", show_default=True)
def build_candidate_cards(config_path: str) -> None:
    summary = run_build_candidate_cards(Path(config_path))
    click.echo(f"Stage 12 candidate cards: {summary['cards']}")


@main.command("render-candidate-dossiers")
@click.option("--config", "config_path", default="configs/stage12_candidate_dossiers.yaml", show_default=True)
def render_candidate_dossiers(config_path: str) -> None:
    summary = run_render_candidate_dossiers(Path(config_path))
    click.echo(f"Stage 12 candidate dossiers: {summary['dossiers']}")


@main.command("render-final-figures")
@click.option("--config", "config_path", default="configs/stage12_candidate_dossiers.yaml", show_default=True)
def render_final_figures(config_path: str) -> None:
    summary = run_render_final_figures(Path(config_path))
    click.echo(f"Stage 12 figures: {summary['figures']}")


@main.command("build-provenance-bundle")
@click.option("--config", "config_path", default="configs/stage12_candidate_dossiers.yaml", show_default=True)
def build_provenance_bundle(config_path: str) -> None:
    summary = run_build_provenance_bundle(Path(config_path))
    click.echo(f"Stage 12 provenance manifest: {summary['manifest']}")


@main.command("build-model-card")
@click.option("--config", "config_path", default="configs/stage12_candidate_dossiers.yaml", show_default=True)
def build_model_card(config_path: str) -> None:
    summary = run_build_model_card(Path(config_path))
    click.echo(f"Stage 12 model card: {summary['model_card']}")


@main.command("build-dataset-card")
@click.option("--config", "config_path", default="configs/stage12_candidate_dossiers.yaml", show_default=True)
def build_dataset_card(config_path: str) -> None:
    summary = run_build_dataset_card(Path(config_path))
    click.echo(f"Stage 12 dataset card: {summary['dataset_card']}")


@main.command("report-stage12")
@click.option("--config", "config_path", default="configs/stage12_candidate_dossiers.yaml", show_default=True)
def report_stage12(config_path: str) -> None:
    summary = run_report_stage12(Path(config_path))
    click.echo(f"Stage 12 report: {summary['report']}")


@main.command("run-stage12")
@click.option("--config", "config_path", default="configs/stage12_candidate_dossiers.yaml", show_default=True)
def run_stage12(config_path: str) -> None:
    summary = run_stage12_all(Path(config_path))
    click.echo(f"Stage 12 complete: {summary['report']}; manifest: {summary['manifest']}")


@main.command("fetch-structures")
@click.option("--config", "project_config", default="configs/project.yaml", show_default=True)
@click.option("--stage-config", default="configs/data_sources.yaml", show_default=True)
@click.option("--out", "out_dir", default="data/processed/stage1", show_default=True)
@click.option("--workers", default=None, type=int, help="Download/extraction workers.")
@click.option("--force", is_flag=True, help="Re-download mmCIF files and overwrite derived files.")
def fetch_structures(
    project_config: str,
    stage_config: str,
    out_dir: str,
    workers: int | None,
    force: bool,
) -> None:
    """Build the Stage 1 EGFR co-crystal benchmark."""
    summary = run_stage1(
        project_config=Path(project_config),
        stage_config=Path(stage_config),
        out_dir=Path(out_dir),
        workers=workers,
        force=force,
    )
    click.echo(f"Stage 1 status: {summary['status']}")
    click.echo(f"Structures table: {summary['structures_table']}")
    click.echo(f"Ligands table: {summary['ligands_table']}")


@main.command("prepare-receptors")
def prepare_receptors() -> None:
    raise click.ClickException("Not implemented yet. Stage 1 writes clean receptor PDB files.")


@main.command("prepare-ligands")
def prepare_ligands() -> None:
    raise click.ClickException("Not implemented yet.")


@main.command("redock")
def redock() -> None:
    raise click.ClickException("Not implemented yet.")


@main.command("crossdock")
def crossdock() -> None:
    raise click.ClickException("Not implemented yet.")


@main.command("rescore-gnina")
@click.option("--config", "config_path", default="configs/stage4_rescoring.yaml", show_default=True)
def rescore_gnina(config_path: str) -> None:
    summary = run_stage4_gnina(Path(config_path))
    click.echo(f"GNINA runs: {summary['runs']}; scores: {summary['scores']}")


@main.command("compute-prolif")
def compute_prolif() -> None:
    raise click.ClickException("Not implemented yet.")


@main.command("train-pose-model")
def train_pose_model() -> None:
    raise click.ClickException("Not implemented yet.")


@main.command("build-candidate-library")
def build_candidate_library() -> None:
    raise click.ClickException("Not implemented yet.")


@main.command("screen-candidates")
def screen_candidates() -> None:
    raise click.ClickException("Not implemented yet.")


@main.command("run-agent-loop")
def run_agent_loop() -> None:
    raise click.ClickException("Not implemented yet.")


@main.command("run-md")
def run_md() -> None:
    raise click.ClickException("Not implemented yet.")


@main.command("make-reports")
def make_reports() -> None:
    raise click.ClickException("Not implemented yet.")


if __name__ == "__main__":
    main()
