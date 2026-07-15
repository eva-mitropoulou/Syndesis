.PHONY: test verify

test:
	pytest -q

verify:
	python scripts/verify_results.py
