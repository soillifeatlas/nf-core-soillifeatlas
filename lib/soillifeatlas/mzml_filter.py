"""Keep one MS2 activation type in an mzML file.

Some Orbitrap methods acquire two fragment scans per precursor, for example an ion-trap CID scan
(low mass accuracy) and an Orbitrap HCD scan (high mass accuracy), as on the LTQ Orbitrap Velos
Pro. MZmine merges every fragment scan of a feature, so the low-accuracy scans would blur the
lipid rules and spectral matches. This keeps every MS1 scan and only the fragment scans of the
chosen activation.

A fragment scan's activation is read from its filter string (``@hcd``, ``@cid``) and, when the
file has no filter strings, from the activation term of its precursor. The output is a plain
(non-indexed) mzML: the index of the input would no longer point at the right bytes.
"""
from __future__ import annotations

import re

_SPECTRUM = re.compile(r'<spectrum [\s\S]*?</spectrum>')
_LEVEL = re.compile(r'value="(\d+)" name="ms level"')
_FILTER = re.compile(r'value="([^"]*)" name="filter string"')
_TERMS = {
    'hcd': 'name="beam-type collision-induced dissociation"',
    'cid': 'name="collision-induced dissociation"',
}


def _activation_matches(block: str, activation: str) -> bool:
    filt = _FILTER.search(block)
    if filt and '@' in filt.group(1):
        return f'@{activation}' in filt.group(1).lower()
    return _TERMS[activation] in block


def keep_activation(text: str, activation: str) -> tuple[str, dict[str, int]]:
    """Return the mzML text with only MS1 scans and fragment scans of ``activation`` ('hcd' or 'cid'),
    and the counts of MS1 scans, fragment scans kept and fragment scans dropped."""
    activation = activation.lower()
    if activation not in _TERMS:
        raise ValueError(f"activation must be one of {sorted(_TERMS)}, not {activation!r}")
    counts = {'ms1': 0, 'kept_msn': 0, 'dropped_msn': 0}
    kept: list[str] = []
    for block in _SPECTRUM.findall(text):
        level = _LEVEL.search(block)
        if level is None or level.group(1) == '1':
            counts['ms1'] += 1
            kept.append(block)
        elif _activation_matches(block, activation):
            counts['kept_msn'] += 1
            kept.append(block)
        else:
            counts['dropped_msn'] += 1
    kept = [re.sub(r'(<spectrum [^>]*?)index="\d+"', rf'\g<1>index="{i}"', b, count=1) for i, b in enumerate(kept)]

    first = text.index('<spectrum ') if kept or '<spectrum ' in text else None
    head_end = text.rfind('<spectrumList', 0, first)
    head = text[:head_end] + re.sub(r'count="\d+"', f'count="{len(kept)}"', text[head_end:first], count=1)
    tail = text[text.index('</spectrumList>'):text.index('</mzML>') + len('</mzML>')]
    # Drop the indexedmzML wrapper: keep the XML declaration and start at <mzML>.
    decl = '<?xml version="1.0" encoding="utf-8"?>\n'
    head = decl + head[head.index('<mzML'):]
    return head + ''.join(kept) + tail + '\n', counts
