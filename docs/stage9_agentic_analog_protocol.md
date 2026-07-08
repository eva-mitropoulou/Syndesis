# Stage 9 Agentic Analog Protocol

Stage 9 performs interaction-constrained analog optimization from Stage 8 seeds. Agents may propose transformation requests, but RDKit creates molecules and the existing Stage 8 mini-screen performs docking, GNINA rescoring, ProLIF interaction checks, and Stage 6 pose-confidence application.

The stage does not claim experimental activity and does not accept free-form LLM molecules as final molecules.

Default run order:

1. Select Stage 8 seeds.
2. Detect protected and editable sites.
3. Write the allowed transformation library.
4. Enumerate deterministic RDKit analogs.
5. Record REINVENT4 and local LLM strategy status.
6. Validate molecules and medchem scope.
7. Run Stage 8 mini-screen for hard-scope analogs.
8. Apply binding-mode-preserving acceptance policy.
9. Benchmark strategies and write the Stage 9 report.
