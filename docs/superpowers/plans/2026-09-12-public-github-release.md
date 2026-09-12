# Public GitHub Release Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make BGC-XPLORER safe to publish, runnable from a versioned GHCR image, and continuously verified and released by GitHub Actions.

**Architecture:** The public Compose file pulls a tagged GHCR image and persists inputs, results, configuration, and databases as host directories. A development override builds the same image locally. CI validates lightweight behavior on pull requests, while release tags build and publish the immutable application image without databases.

**Tech Stack:** Docker/OCI, Docker Compose, GitHub Actions, GHCR, Python unittest, Snakemake mock mode

**Spec:** Approved deployment design in the project conversation, 2026-09-12

## Global Constraints

- NVIDIA credentials must never be committed or printed by CI.
- Reference databases, analysis inputs, results, caches, and local environments remain outside Git and outside the image.
- Linux amd64 is the first supported release platform.
- Tool versions and the ARTS source revision are immutable for a tagged release.
- The public startup command is `docker compose up -d` after creating `.env`.

---

### Task 1: Define and verify the public source boundary

**Files:**
- Modify: `.gitignore`
- Modify: `.dockerignore`
- Create: `scripts/check_public_tree.py`
- Create: `tests/test_check_public_tree.py`

**Interfaces:**
- Consumes: repository-relative file paths
- Produces: `find_public_tree_violations(paths) -> list[str]`

- [ ] Write tests proving secrets, databases, results, nested Git metadata, and runtime caches are rejected.
- [ ] Run `python3 -m unittest tests.test_check_public_tree -v` and confirm the missing module failure.
- [ ] Implement the public-tree validator and ignore rules.
- [ ] Run the test and scan the candidate public file list.

### Task 2: Separate production deployment from local development

**Files:**
- Modify: `docker-compose.yml`
- Create: `docker-compose.dev.yml`
- Modify: `.env.example`
- Modify: `README.md`
- Modify: `Dockerfile`

**Interfaces:**
- Consumes: `BGC_IMAGE`, `BGC_VERSION`, `BGC_PORT`, `NVIDIA_API_KEY`, `NVIDIA_MODEL`
- Produces: production pull deployment and explicit local build deployment

- [ ] Set the production service image to `ghcr.io/${BGC_IMAGE}:${BGC_VERSION}` without `build:`.
- [ ] Add a development override containing `build: .` and the local image tag.
- [ ] Pin application and analysis tool versions used by the measured working image.
- [ ] Exclude ARTS Git metadata, bundled references, uploads, and databases from the build context.
- [ ] Validate both Compose configurations with a placeholder API key.

### Task 3: Add CI and tagged GHCR publishing

**Files:**
- Create: `.github/workflows/ci.yml`
- Create: `.github/workflows/release-container.yml`
- Create: `tests/test_release_metadata.py`

**Interfaces:**
- Consumes: pushes, pull requests, and tags matching `v*.*.*`
- Produces: test results and linux/amd64 images tagged with semantic version, major/minor aliases, and commit SHA

- [ ] Write metadata tests that require least-privilege permissions, pinned action majors, linux/amd64, and GHCR tags.
- [ ] Run the metadata test and confirm it fails while workflows are absent.
- [ ] Add CI checks for unit tests, Python compilation, shell syntax, public-tree validation, Compose validation, and Snakemake mock dry-run.
- [ ] Add the release workflow using GitHub OIDC-provided `GITHUB_TOKEN` with `packages: write`.
- [ ] Run all local lightweight checks.

### Task 4: Initialize the publication branch and review the deliverable

**Files:**
- Create: `.git/` metadata only

**Interfaces:**
- Consumes: cleaned public source tree
- Produces: Git branch `release/public-github` ready for remote configuration and review

- [ ] Initialize Git without adding a remote and create `release/public-github`.
- [ ] Inspect `git status --short --ignored` and verify large/generated directories are ignored.
- [ ] Run the complete verification suite and report any check that requires GitHub infrastructure.
