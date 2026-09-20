# Shared helpers for talking to the GitHub release API without the gh CLI.
# Sourced by fetch_sdk.sh and upload_sdk.sh; not executable on its own.

# The manylinux build image ships no `python3` on PATH, only the per-tag
# interpreters under /opt/python. Find whichever exists.
msani_python() {
    if command -v python3 >/dev/null 2>&1; then
        command -v python3
        return
    fi
    local candidate
    for candidate in /opt/python/cp312-cp312/bin/python3 /opt/python/cp31*/bin/python3 \
                     /usr/libexec/platform-python; do
        if [ -x "$candidate" ]; then
            printf '%s\n' "$candidate"
            return
        fi
    done
    printf 'No python3 available for JSON parsing.\n' >&2
    return 1
}

msani_api() {
    curl -sfL -H "Authorization: Bearer $token" \
        -H 'Accept: application/vnd.github+json' \
        -H 'X-GitHub-Api-Version: 2022-11-28' "$@"
}
