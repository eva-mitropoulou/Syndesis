# Stage 9 Local LLM Setup

The default configuration expects a local Qwen3-32B quantized model served through the configured local provider. This workspace is configured for `Qwen/Qwen3-32B-AWQ` through vLLM at `http://localhost:8000`.

Supported providers are:

- `ollama`
- `llama_cpp_server`
- `vllm`
- `openai_compatible_local`
- `disabled_mock_for_tests`

Production runs must use a real local provider. If configured LLM strategies are enabled and the provider is unavailable, Stage 9 fails instead of recording fallback proposals.
