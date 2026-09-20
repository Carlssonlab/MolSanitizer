#!/usr/bin/env bash
# Attach one SDK tarball to the SDK pre-release, replacing any existing asset
# of the same name.
#
# Mirrors fetch_sdk.sh: uses the gh CLI when present, and the REST API
# otherwise, because the manylinux build container has no gh package.
#
# Usage: upload_sdk.sh <release-tag> <file>
# Requires GH_TOKEN (or GITHUB_TOKEN) with write access to the repository.
set -euo pipefail

tag=${1:?Usage: upload_sdk.sh <release-tag> <file>}
file=${2:?Usage: upload_sdk.sh <release-tag> <file>}
repo=${GITHUB_REPOSITORY:?GITHUB_REPOSITORY must be set}
api=${GITHUB_API_URL:-https://api.github.com}
uploads=${GITHUB_SERVER_URL:+https://uploads.github.com}
uploads=${uploads:-https://uploads.github.com}
token=${GH_TOKEN:-${GITHUB_TOKEN:-}}
name=$(basename "$file")
test -s "$file"

# shellcheck source=/dev/null
. "$(dirname "${BASH_SOURCE[0]}")/_gh_api.sh"

if command -v gh >/dev/null 2>&1; then
    gh release upload "$tag" "$file" --repo "$repo" --clobber
else
    : "${token:?GH_TOKEN or GITHUB_TOKEN is required without the gh CLI}"
    python=$(msani_python)
    release=$(msani_api "$api/repos/$repo/releases/tags/$tag")
    read -r release_id asset_id <<EOF
$(printf '%s' "$release" | "$python" -c '
import json, sys
name = sys.argv[1]
release = json.load(sys.stdin)
existing = [a["id"] for a in release["assets"] if a["name"] == name]
print(release["id"], existing[0] if existing else "")
' "$name")
EOF
    if [ -n "${asset_id:-}" ]; then
        printf 'Replacing existing asset %s\n' "$name"
        curl -sfL -X DELETE -H "Authorization: Bearer $token" \
            -H 'X-GitHub-Api-Version: 2022-11-28' \
            "$api/repos/$repo/releases/assets/$asset_id"
    fi
    curl -fL --retry 3 -X POST \
        -H "Authorization: Bearer $token" \
        -H 'Content-Type: application/gzip' \
        -H 'X-GitHub-Api-Version: 2022-11-28' \
        --data-binary "@$file" \
        "$uploads/repos/$repo/releases/$release_id/assets?name=$name" >/dev/null
fi

printf 'Uploaded %s to %s\n' "$name" "$tag"
