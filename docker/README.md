# CapScreen Docker images

CapScreen uses two images:

| Image | Registry tag | Contents | Rebuild when |
| --- | --- | --- | --- |
| Base | `ghcr.io/hassansaei/capscreen-base:latest` (no version) | OS packages, mamba, fastp, conda envs | `conda/*.yml`, system tools, or `docker/DockerFile.base` change |
| App (final) | `saei/capscreen:<version>` and optionally `:latest` | CapScreen source + CLI | Python/scripts/package change |

The base image is unversioned (`:latest` only) and published to GitHub Packages (GHCR).
The final app image carries the CapScreen version tag.

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

## When to rebuild the final (app) image

Rebuild the app image when CapScreen code changes (scripts, package, templates),
or after a new base was published and you want the app to pick it up.

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
