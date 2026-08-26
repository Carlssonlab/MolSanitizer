"""Utilities for transparently detecting compressed inputs."""


_COMPRESSION_SIGNATURES = (
    (b'\x1f\x8b', 'gzip'),
    (b'\xfd7zXZ\x00', 'xz'),
)


def detect_input_compression(input_file):
    """Return the compression type detected from a file signature."""
    with open(input_file, 'rb') as stream:
        signature = stream.read(max(len(magic) for magic, _ in _COMPRESSION_SIGNATURES))

    for magic, compression in _COMPRESSION_SIGNATURES:
        if signature.startswith(magic):
            return compression
    return None
