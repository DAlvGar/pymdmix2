#!/usr/bin/env bash
# Build the pyMDmix2 Docker image.
#
# Usage:
#   ./scripts/build_docker.sh [OPTIONS]
#
# Environment variables (all optional):
#   IMAGE_NAME   Docker image name  (default: pymdmix)
#   IMAGE_TAG    Docker image tag   (default: latest)
#   NO_CACHE     Set to 1 to pass --no-cache to docker build
#
# Examples:
#   ./scripts/build_docker.sh
#   IMAGE_TAG=0.3.0 ./scripts/build_docker.sh
#   NO_CACHE=1 IMAGE_TAG=dev ./scripts/build_docker.sh

set -euo pipefail

IMAGE_NAME="${IMAGE_NAME:-pymdmix}"
IMAGE_TAG="${IMAGE_TAG:-latest}"
NO_CACHE="${NO_CACHE:-0}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

EXTRA_ARGS=()
if [[ "${NO_CACHE}" == "1" ]]; then
    EXTRA_ARGS+=("--no-cache")
fi

echo "==> Building Docker image: ${IMAGE_NAME}:${IMAGE_TAG}"
echo "    Context: ${REPO_ROOT}"
echo "    Dockerfile: ${REPO_ROOT}/Dockerfile"
echo ""

docker build \
    --tag "${IMAGE_NAME}:${IMAGE_TAG}" \
    --file "${REPO_ROOT}/Dockerfile" \
    "${EXTRA_ARGS[@]}" \
    "${REPO_ROOT}"

echo ""
echo "==> Build complete: ${IMAGE_NAME}:${IMAGE_TAG}"
echo ""
echo "    Run with:"
echo "      docker run --rm -v \$(pwd):/work ${IMAGE_NAME}:${IMAGE_TAG} --help"
echo ""
echo "    Push to registry:"
echo "      docker push ${IMAGE_NAME}:${IMAGE_TAG}"
