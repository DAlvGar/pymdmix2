#!/usr/bin/env bash
# Build a Singularity / Apptainer image for pyMDmix2.
#
# The script tries three approaches in order:
#   1. If a local Docker image exists, convert it with docker-daemon://
#   2. If a Singularity definition file (Singularity.def) is present, build from it.
#   3. Pull directly from a remote Docker registry (requires IMAGE_REGISTRY to be set).
#
# Usage:
#   ./scripts/build_singularity.sh [OPTIONS]
#
# Environment variables (all optional):
#   IMAGE_NAME      Docker image name used as source  (default: pymdmix)
#   IMAGE_TAG       Docker image tag                  (default: latest)
#   SIF_FILE        Output .sif filename              (default: pymdmix.sif)
#   IMAGE_REGISTRY  Remote registry, e.g. docker://ghcr.io/org/pymdmix
#   DEF_FILE        Path to a .def file to build from (overrides auto-detect)
#   FAKEROOT        Set to 1 to pass --fakeroot (needed on systems without root)
#
# Examples:
#   ./scripts/build_singularity.sh
#   SIF_FILE=pymdmix-0.3.0.sif ./scripts/build_singularity.sh
#   FAKEROOT=1 DEF_FILE=Singularity.def ./scripts/build_singularity.sh
#   IMAGE_REGISTRY=docker://ghcr.io/DAlvGar/pymdmix ./scripts/build_singularity.sh

set -euo pipefail

IMAGE_NAME="${IMAGE_NAME:-pymdmix}"
IMAGE_TAG="${IMAGE_TAG:-latest}"
SIF_FILE="${SIF_FILE:-pymdmix.sif}"
IMAGE_REGISTRY="${IMAGE_REGISTRY:-}"
DEF_FILE="${DEF_FILE:-}"
FAKEROOT="${FAKEROOT:-0}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# ---------------------------------------------------------------------------
# Detect singularity or apptainer
# ---------------------------------------------------------------------------
if command -v apptainer &>/dev/null; then
    SING=apptainer
elif command -v singularity &>/dev/null; then
    SING=singularity
else
    echo "ERROR: Neither 'apptainer' nor 'singularity' found on PATH." >&2
    echo "       Install one of them: https://apptainer.org/docs/admin/main/installation.html" >&2
    exit 1
fi

EXTRA_ARGS=()
if [[ "${FAKEROOT}" == "1" ]]; then
    EXTRA_ARGS+=("--fakeroot")
fi

# ---------------------------------------------------------------------------
# Choose build source
# ---------------------------------------------------------------------------
if [[ -n "${DEF_FILE}" ]]; then
    # Explicit definition file provided
    echo "==> Building from definition file: ${DEF_FILE}"
    "${SING}" build "${EXTRA_ARGS[@]}" "${SIF_FILE}" "${DEF_FILE}"

elif docker image inspect "${IMAGE_NAME}:${IMAGE_TAG}" &>/dev/null; then
    # Local Docker image found — convert directly (no daemon push needed)
    echo "==> Converting local Docker image ${IMAGE_NAME}:${IMAGE_TAG} → ${SIF_FILE}"
    "${SING}" build "${EXTRA_ARGS[@]}" "${SIF_FILE}" "docker-daemon://${IMAGE_NAME}:${IMAGE_TAG}"

elif [[ -f "${REPO_ROOT}/Singularity.def" ]]; then
    # Fall back to the repo's definition file
    echo "==> Building from ${REPO_ROOT}/Singularity.def"
    "${SING}" build "${EXTRA_ARGS[@]}" "${SIF_FILE}" "${REPO_ROOT}/Singularity.def"

elif [[ -n "${IMAGE_REGISTRY}" ]]; then
    # Pull from a remote registry
    echo "==> Pulling from remote registry: ${IMAGE_REGISTRY}"
    "${SING}" build "${EXTRA_ARGS[@]}" "${SIF_FILE}" "${IMAGE_REGISTRY}"

else
    echo "ERROR: No build source found. Provide one of:" >&2
    echo "       - A local Docker image '${IMAGE_NAME}:${IMAGE_TAG}' (run build_docker.sh first)" >&2
    echo "       - DEF_FILE=<path/to/definition.def>" >&2
    echo "       - IMAGE_REGISTRY=docker://<registry/image:tag>" >&2
    echo "       - A Singularity.def file at the repository root" >&2
    exit 1
fi

echo ""
echo "==> Singularity image written to: ${SIF_FILE}"
echo ""
echo "    Run with:"
echo "      ${SING} run ${SIF_FILE} --help"
echo "      ${SING} run --bind \$(pwd):/work --pwd /work ${SIF_FILE} [COMMAND]"
