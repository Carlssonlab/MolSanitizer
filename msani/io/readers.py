"""Utilities for detecting input compression and column separators."""

import bz2
import gzip
import lzma


_COMPRESSION_SIGNATURES = (
    (b'\x1f\x8b', 'gzip'),
    (b'BZh', 'bz2'),
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


def detect_input_separator(input_file, compression):
    """Inspect at most 64 KiB of decompressed input, leaving parsing to pandas.

    The first nonblank record determines the separator for the entire file.
    Tabs preserve optional CXSMILES extensions even if they occur much later.
    """
    opener = {'gzip': gzip.open, 'bz2': bz2.open, 'xz': lzma.open}.get(compression, open)
    sample_size = 64 * 1024
    with opener(input_file, 'rb') as stream:
        sample = stream.read(sample_size)
    lines = sample.split(b'\n')
    if len(sample) == sample_size:
        # Do not infer a delimiter from a potentially incomplete final record.
        lines.pop()
    for line in lines:
        if line.strip():
            # A tab after an existing space-separated ID is an extra column,
            # not the structure/ID delimiter.
            first_field, tab, _ = line.partition(b'\t')
            fields = first_field.strip().split(maxsplit=1)
            if tab and fields and (len(fields) == 1 or fields[1].startswith(b'|')):
                return '\t'
            return r'\s+'
    if len(sample) == sample_size:
        raise ValueError(
            f'{input_file}: no complete input record within the first 64 KiB. '
            'Remove leading blank lines or use -e for tab-separated input.'
        )
    return r'\s+'


def validate_input_chunk(chunk, input_file, record_offset=0):
    """Catch delimiter mistakes before a chunk reaches molecular processing."""
    ids = chunk['ids']
    invalid = ids.str.match(r'\s*(?:\||$)', na=True)
    if invalid.any():
        position = invalid.to_numpy().nonzero()[0][0]
        value = ids.iloc[position]
        reason = ('CXSMILES extension in the ID column'
                  if isinstance(value, str) and value.lstrip().startswith('|')
                  else 'missing compound ID')
        raise ValueError(
            f'{input_file}: input record {record_offset + position + 1}: {reason}. '
            'Use a consistent separator throughout the file. For CXSMILES, '
            'use tabs between structures and IDs, or double-quote the complete '
            'CXSMILES field for space-separated input and omit -e.'
        )
