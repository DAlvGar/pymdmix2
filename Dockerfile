# syntax=docker/dockerfile:1
#
# pyMDmix2 Docker image
#
# Provides AmberTools (cpptraj, LEaP, sander/pmemd via host or GPU passthrough)
# plus the full pymdmix Python package with all optional dependencies.
#
# Build:
#   docker build -t pymdmix:latest .
#
# Run:
#   docker run --rm -v $(pwd):/work pymdmix:latest --help

FROM continuumio/miniconda3:24.1.2-0

LABEL org.opencontainers.image.title="pyMDmix2" \
      org.opencontainers.image.description="Molecular Dynamics with organic solvent mixtures for drug discovery" \
      org.opencontainers.image.authors="Daniel Alvarez-Garcia <algarcia.daniel@gmail.com>" \
      org.opencontainers.image.source="https://github.com/DAlvGar/pymdmix2" \
      org.opencontainers.image.licenses="MIT"

# System build dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        git && \
    rm -rf /var/lib/apt/lists/*

# AmberTools + Python 3.11 via conda-forge.
# ambertools on conda-forge includes cpptraj, LEaP, sander, antechamber, etc.
RUN conda install -y -c conda-forge \
        python=3.11 \
        ambertools && \
    conda clean -afy

# Expose AMBERHOME so pymdmix can locate tleap, cpptraj, and force-field data
# when running in non-interactive shells (Docker entrypoint, CI runners, etc.).
ENV AMBERHOME=/opt/conda

# Copy project source and install pymdmix with all optional extras
WORKDIR /opt/pymdmix
COPY . .
RUN pip install --no-cache-dir ".[full]"

# User working directory — mount project data here
WORKDIR /work

ENTRYPOINT ["pymdmix"]
CMD ["--help"]
