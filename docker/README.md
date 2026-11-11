# CapScreen Docker images

CapScreen uses two images:

| Image | Registry tag | Contents | Rebuild when |
| --- | --- | --- | --- |
| Base | `ghcr.io/hassansaei/capscreen-base:latest` (no version) | OS packages, mamba, fastp, conda envs | `conda/*.yml`, system tools, or `docker/DockerFile.base` change |
| App (final) | `ghcr.io/hassansaei/capscreen:<version>` (+ `:latest`); also `saei/capscreen` if Docker Hub secrets are set | CapScreen source + CLI | New GitHub Release (automatic) or code change |

The base image is unversioned (`:latest` only) and published to GitHub Packages (GHCR).
The final app image carries the CapScreen version tag.
CI builds it from the published base via `.github/workflows/docker-app.yml`.

Ownership is limited to input/output/reference dirs in the base image and to
`/opt/capscreen` in the app image. Neither image runs `chown -R /opt`.

## When to rebuild the base

Rebuild and push the base only if you:

- change `conda/capscreen.yml` or `conda/capscreen-stat.yml`
- add/change system packages or tools in `docker/DockerFile.base`
- need a fresh conda solve

Do **not** rebuild the base for script-only edits (`stat.py`, CLI, config, docs).

### Local base build

```bash
docker build --platform linux/amd64 \
  -t ghcr.io/hassansaei/capscreen-base:latest \
  -f docker/DockerFile.base .
```

### Publish base to GitHub Packages

```bash
echo "$GITHUB_TOKEN" | docker login ghcr.io -u USERNAME --password-stdin
docker push ghcr.io/hassansaei/capscreen-base:latest
```

CI publishes the base automatically on changes to `conda/**` or
`docker/DockerFile.base` via `.github/workflows/docker-base.yml`
(also runnable with **Actions → Publish CapScreen base image → Run workflow**).

## Automatic app image on release

Publishing a GitHub Release (including pre-releases such as `2.0.0-alpha.2`)
triggers **Actions → Publish CapScreen app image**. That workflow:

1. Uses `ghcr.io/hassansaei/capscreen-base:latest` as `BASE_IMAGE`
2. Builds `docker/DockerFile` for `linux/amd64`
3. Pushes `ghcr.io/hassansaei/capscreen:<tag>` and `:latest`
4. Converts the pushed image to Apptainer/Singularity `capscreen-<tag>.sif`
5. Attaches that `.sif` to the GitHub Release (and uploads it as a workflow artifact)

Optional Docker Hub publish (`saei/capscreen:<tag>` and `:latest`) runs when
these repository secrets exist:

- `DOCKERHUB_USERNAME`
- `DOCKERHUB_TOKEN` (access token, not password)

Add them under **Settings → Secrets and variables → Actions**.

You can also run the workflow manually:
**Actions → Publish CapScreen app image → Run workflow**
(optionally pass a version; otherwise it reads `capscreen/version.py`).
Manual runs still build the `.sif` and upload it as an Actions artifact.

HPC / Apptainer example (after downloading the release asset):

```bash
apptainer run capscreen-2.2.0.sif --version
apptainer run capscreen-2.2.0.sif pipeline --vp-seq ...
```

Release checklist:

1. Push commits that bump `capscreen/version.py` (post-commit hook for alphas)
2. Create a GitHub Release whose tag matches that version (e.g. `2.0.0-alpha.2`)
3. Wait for the app image workflow to finish; pull from GHCR (and Docker Hub if configured)
4. Download `capscreen-<tag>.sif` from the release assets if you need Singularity/Apptainer

## When to rebuild the final (app) image locally

Rebuild locally when iterating on CapScreen code before a release,
or after a new base was published and you want to test against it.

### Local app build (uses published base)

```bash
VERSION=$(python -c "import re; t=open('capscreen/version.py').read(); print(re.search(r\"__version__\\s*=\\s*['\\\"]([^'\\\"]+)['\\\"]\", t).group(1))")

docker build --platform linux/amd64 \
  -t saei/capscreen:${VERSION} \
  -t saei/capscreen:latest \
  -f docker/DockerFile .
```

### Build app against a local base (before GHCR publish)

```bash
docker build --platform linux/amd64 \
  --build-arg BASE_IMAGE=ghcr.io/hassansaei/capscreen-base:latest \
  -t saei/capscreen:latest \
  -f docker/DockerFile .
```

If the base exists only locally under that same tag, Docker uses it without pulling.

### Push final image to Docker Hub

```bash
docker push saei/capscreen:latest
# optional versioned tag:
# docker push saei/capscreen:${VERSION}
```

## Run

```bash
docker run --rm --platform linux/amd64 \
  -v /path/to/input:/opt/capscreen/input \
  -v /path/to/output:/opt/capscreen/output \
  saei/capscreen:latest \
  --version
```
