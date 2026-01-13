# RNAstructure docker wrappers

These scripts provide drop-in replacements for RNAstructure executables (e.g. `AllSub`, `Partition`) by running them inside a Docker image.

- Configure the image with `RNASTRUCTURE_DOCKER_IMAGE` (default: `rnanneal/rnastructure:v6.4.0`, but the benchmark runner config uses ECR URIs).
- The wrappers translate host file paths into container paths by bind-mounting the relevant parent directories.

Example (manual):

```bash
export RNASTRUCTURE_DOCKER_IMAGE=rnanneal/rnastructure:v6.4.0
benchmark_runner/tools/rnastructure_docker/AllSub in.fa out.ct -a 20
```
