# Stage 9 Agent Roles

Agents are structured prompt templates, not autonomous programs.

- Medicinal Chemist Agent: proposes conservative scaffold-preserving transformations.
- Docking Scientist Agent: checks whether proposed edits are compatible with receptor-state and pose context.
- Interaction Analyst Agent: protects key ProLIF contacts and binding-mode evidence.
- ADMET/Safety Critic Agent: flags medchem risks without clinical claims.
- Skeptical Reviewer Agent: rejects score-hacking and weak evidence.
- Orchestrator: consolidates tool-executable transformation requests and never invents final molecules directly.

Agent output must be JSON only and must contain concise rationale summaries. Chain-of-thought fields are not required and are rejected by schema tests.
