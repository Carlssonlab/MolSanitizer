#!/usr/bin/env bash
# Download one SDK asset from the SDK pre-release into the current directory.
#
# Uses the gh CLI when present (macOS and Windows runners) and falls back to
# the REST API with curl, because the manylinux build container ships neither
# gh nor a package for it.
#
# Usage: fetch_sdk.sh <release-tag> <asset-name>
# Requires GH_TOKEN (or GITHUB_TOKEN) with read access to the repository.
set -euo pipefail

tag=${1:?Usage: fetch_sdk.sh <release-tag> <asset-name>}
asset=${2:?Usage: fetch_sdk.sh <release-tag> <asset-name>}
repo=${GITHUB_REPOSITORY:?GITHUB_REPOSITORY must be set}
api=${GITHUB_API_URL:-https://api.github.com}
token=${GH_TOKEN:-${GITHUB_TOKEN:-}}

# shellcheck source=/dev/null
. "$(dirname "${BASH_SOURCE[0]}")/_gh_api.sh"

if command -v gh >/dev/null 2>&1; then
    gh release download "$tag" --repo "$repo" --pattern "$asset" --clobber
else
    : "${token:?GH_TOKEN or GITHUB_TOKEN is required without the gh CLI}"
    python=$(msani_python)
    url=$(msani_api "$api/repos/$repo/releases/tags/$tag" | "$python" -c '
import json, sys
name = sys.argv[1]
assets = json.load(sys.stdin)["assets"]
matches = [a["url"] for a in assets if a["name"] == name]
if len(matches) != 1:
    sys.exit("Expected one asset named %s; release holds: %s"
             % (name, ", ".join(a["name"] for a in assets) or "nothing"))
print(matches[0])
' "$asset")
    # The asset API returns the file itself only with an octet-stream Accept.
    curl -sfL --retry 3 -H "Authorization: Bearer $token" \
        -H 'Accept: application/octet-stream' "$url" -o "$asset"
fi

test -s "$asset"
printf 'Fetched %s (%s bytes)\n' "$asset" "$(wc -c < "$asset")"
