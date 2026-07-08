.PHONY: stage1 stage2 stage3 test lint

stage1:
	egfrforge build-cocrystal-benchmark --config configs/stage1_cocrystal_benchmark.yaml

stage2:
	egfrforge build-receptor-features --config configs/stage2_receptor_ensemble.yaml
	egfrforge cluster-receptors --config configs/stage2_receptor_ensemble.yaml
	egfrforge select-receptor-ensemble --config configs/stage2_receptor_ensemble.yaml
	egfrforge report-stage2 --config configs/stage2_receptor_ensemble.yaml

stage3:
	egfrforge prepare-docking-inputs --config configs/stage3_redocking_crossdocking.yaml
	egfrforge build-docking-task-matrix --config configs/stage3_redocking_crossdocking.yaml
	egfrforge run-redocking --config configs/stage3_redocking_crossdocking.yaml
	egfrforge compute-docking-rmsd --config configs/stage3_redocking_crossdocking.yaml
	egfrforge run-pose-sanity-checks --config configs/stage3_redocking_crossdocking.yaml
	egfrforge label-stage3-poses --config configs/stage3_redocking_crossdocking.yaml
	egfrforge report-stage3 --config configs/stage3_redocking_crossdocking.yaml

test:
	pytest -q

lint:
	ruff check src tests
