#!/usr/bin/env bash
# Run AmberTools and real-data end-to-end tests inside the pyMDmix2 Docker image.
#
# The script builds (or reuses) the Docker image and executes the tests
# tagged with @pytest.mark.ambertools, @pytest.mark.cpptraj, and/or
# @pytest.mark.real_data so the regular (non-AmberTools) test suite is
# not re-run.
#
# Test tiers
# ----------
# ambertools   – requires tleap / LEaP (no real trajectory data needed)
# cpptraj      – requires cpptraj binary (subset of ambertools tests)
# real_data    – requires the Git LFS binary files in tests/data/
#                (traj.nc, pep_WAT_WAT_1.pdb, ETA_CT.dx)
#
# The MARKER variable (default "ambertools") selects which tests to run.
# Use a pytest marker *expression* to run multiple tiers at once, e.g.:
#   MARKER="ambertools or real_data"   – all AmberTools + all real-data tests
#   MARKER="cpptraj and real_data"     – cpptraj tests on the real trajectory
#
# Usage:
#   ./scripts/run_ambertools_tests.sh [OPTIONS]
#
# Environment variables (all optional):
#   IMAGE_NAME    Docker image name   (default: pymdmix)
#   IMAGE_TAG     Docker image tag    (default: latest)
#   REBUILD       Set to 1 to force a fresh image build (default: 0)
#   PYTEST_ARGS   Extra args forwarded to pytest (default: "-v --tb=short")
#   MARKER        pytest marker expression  (default: "ambertools")
#
# Examples:
#   # Basic run — AmberTools tests only
#   ./scripts/run_ambertools_tests.sh
#
#   # Run all AmberTools + real-data tests together
#   MARKER="ambertools or real_data" ./scripts/run_ambertools_tests.sh
#
#   # Run only the cpptraj tests that use the real 20-frame trajectory
#   MARKER="cpptraj and real_data" ./scripts/run_ambertools_tests.sh
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
# Warn if tests/data/ LFS objects are not checked out
# ---------------------------------------------------------------------------
TRAJ_NC="${REPO_ROOT}/tests/data/amber/traj.nc"
if [[ -f "${TRAJ_NC}" ]]; then
    TRAJ_SIZE=$(wc -c < "${TRAJ_NC}")
    if [[ "${TRAJ_SIZE}" -lt 100000 ]]; then
        echo "WARNING: tests/data/amber/traj.nc looks like an LFS pointer (${TRAJ_SIZE} bytes)."
        echo "         Real-data tests will be skipped inside the container unless you first run:"
        echo "           git lfs pull"
        echo ""
    fi
else
    echo "WARNING: tests/data/amber/traj.nc not found."
    echo "         Real-data tests will be skipped inside the container."
    echo ""
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
echo "==> Running tests (marker: ${MARKER})"
echo "    pytest args: ${PYTEST_ARGS}"
echo ""

# Mount the local tests/ directory so edits are reflected without a rebuild
# and so the real data files in tests/data/ (Git LFS) are available inside
# the container.  The package itself comes from the installed copy inside the
# image.
docker run --rm \
    --name "pymdmix_e2e_$$" \
    -v "${REPO_ROOT}/tests:/opt/pymdmix/tests:ro" \
    --entrypoint pytest \
    "${IMAGE_NAME}:${IMAGE_TAG}" \
    -m "${MARKER}" \
    ${PYTEST_ARGS} \
    /opt/pymdmix/tests/

echo ""
echo "==> Tests finished."
