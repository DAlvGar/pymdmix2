#!/usr/bin/env bash
# Run the AmberTools end-to-end tests inside the pyMDmix2 Docker image.
#
# The script builds (or reuses) the Docker image and executes only the tests
# tagged @pytest.mark.ambertools / @pytest.mark.cpptraj so the regular
# (non-AmberTools) test suite is not re-run.
#
# Usage:
#   ./scripts/run_ambertools_tests.sh [OPTIONS]
#
# Environment variables (all optional):
#   IMAGE_NAME    Docker image name   (default: pymdmix)
#   IMAGE_TAG     Docker image tag    (default: latest)
#   REBUILD       Set to 1 to force a fresh image build (default: 0)
#   PYTEST_ARGS   Extra args forwarded to pytest
#                 (default: "-v --tb=short")
#   MARKER        pytest marker expression  (default: "ambertools")
#
# Examples:
#   # Basic run
#   ./scripts/run_ambertools_tests.sh
#
#   # Run only cpptraj tests, with verbose output
#   MARKER=cpptraj PYTEST_ARGS="-v" ./scripts/run_ambertools_tests.sh
#
#   # Force image rebuild
#   REBUILD=1 ./scripts/run_ambertools_tests.sh
#
#   # Run a single test class
#   PYTEST_ARGS="-v -k TestLeapSession" ./scripts/run_ambertools_tests.sh

set -euo pipefail

IMAGE_NAME="${IMAGE_NAME:-pymdmix}"
IMAGE_TAG="${IMAGE_TAG:-latest}"
REBUILD="${REBUILD:-0}"
PYTEST_ARGS="${PYTEST_ARGS:--v --tb=short}"
MARKER="${MARKER:-ambertools}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# ---------------------------------------------------------------------------
# Check that Docker is available
# ---------------------------------------------------------------------------
if ! command -v docker &>/dev/null; then
    echo "ERROR: 'docker' not found on PATH." >&2
    echo "       Install Docker: https://docs.docker.com/get-docker/" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Build the image if it doesn't exist yet, or if REBUILD=1
# ---------------------------------------------------------------------------
IMAGE_EXISTS=0
if docker image inspect "${IMAGE_NAME}:${IMAGE_TAG}" &>/dev/null; then
    IMAGE_EXISTS=1
fi

if [[ "${REBUILD}" == "1" || "${IMAGE_EXISTS}" == "0" ]]; then
    echo "==> Building Docker image: ${IMAGE_NAME}:${IMAGE_TAG}"
    docker build \
        --tag "${IMAGE_NAME}:${IMAGE_TAG}" \
        --file "${REPO_ROOT}/Dockerfile" \
        "${REPO_ROOT}"
else
    echo "==> Using existing Docker image: ${IMAGE_NAME}:${IMAGE_TAG}"
    echo "    (Set REBUILD=1 to force a rebuild)"
fi

# ---------------------------------------------------------------------------
# Run the tests
# ---------------------------------------------------------------------------
echo ""
echo "==> Running AmberTools e2e tests (marker: ${MARKER})"
echo "    pytest args: ${PYTEST_ARGS}"
echo ""

# Mount the local tests/ directory so edits are reflected without a rebuild.
# The package itself comes from the installed copy inside the image.
docker run --rm \
    --name "pymdmix_e2e_$$" \
    -v "${REPO_ROOT}/tests:/opt/pymdmix/tests:ro" \
    --entrypoint pytest \
    "${IMAGE_NAME}:${IMAGE_TAG}" \
    -m "${MARKER}" \
    ${PYTEST_ARGS} \
    /opt/pymdmix/tests/test_ambertools_e2e.py

echo ""
echo "==> Tests finished."
