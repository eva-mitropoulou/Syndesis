.PHONY: stage1 stage2 stage3 test lint

stage1:
	syndesis build-cocrystal-benchmark --config configs/stage1_cocrystal_benchmark.yaml

stage2:
	syndesis build-receptor-features --config configs/stage2_receptor_ensemble.yaml
	syndesis cluster-receptors --config configs/stage2_receptor_ensemble.yaml
	syndesis select-receptor-ensemble --config configs/stage2_receptor_ensemble.yaml
	syndesis report-stage2 --config configs/stage2_receptor_ensemble.yaml

stage3:
	syndesis prepare-docking-inputs --config configs/stage3_redocking_crossdocking.yaml
	syndesis build-docking-task-matrix --config configs/stage3_redocking_crossdocking.yaml
	syndesis run-redocking --config configs/stage3_redocking_crossdocking.yaml
	syndesis compute-docking-rmsd --config configs/stage3_redocking_crossdocking.yaml
	syndesis run-pose-sanity-checks --config configs/stage3_redocking_crossdocking.yaml
	syndesis label-stage3-poses --config configs/stage3_redocking_crossdocking.yaml
	syndesis report-stage3 --config configs/stage3_redocking_crossdocking.yaml

test:
	pytest -q

lint:
	ruff check src tests
