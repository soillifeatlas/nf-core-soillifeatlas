"""Unit tests for keeping one MS2 activation type in an mzML file (soillifeatlas.mzml_filter)."""
from __future__ import annotations

import re

from soillifeatlas.mzml_filter import keep_activation


def spectrum(scan: int, level: int, filt: str | None, activation: str | None) -> str:
    parts = [f'<spectrum id="controllerType=0 controllerNumber=1 scan={scan}" index="{scan - 1}" defaultArrayLength="1">',
             f'<cvParam cvRef="MS" accession="MS:1000511" value="{level}" name="ms level" />',
             '<scanList count="1"><scan>']
    if filt:
        parts.append(f'<cvParam cvRef="MS" accession="MS:1000512" value="{filt}" name="filter string" />')
    parts.append('</scan></scanList>')
    if activation:
        parts.append(f'<precursorList count="1"><precursor><activation><cvParam cvRef="MS" accession="X" value="" name="{activation}" /></activation></precursor></precursorList>')
    parts.append('<binaryDataArrayList count="0"></binaryDataArrayList></spectrum>')
    return ''.join(parts)


def mzml(spectra: list[str]) -> str:
    return ('<?xml version="1.0" encoding="utf-8"?>\n<indexedmzML xmlns="http://psi.hupo.org/ms/mzml">\n'
            '<mzML xmlns="http://psi.hupo.org/ms/mzml" version="1.1.0"><run id="r" startTimeStamp="2020-02-19T03:14:17Z">'
            f'<spectrumList count="{len(spectra)}" defaultDataProcessingRef="p">' + ''.join(spectra) +
            '</spectrumList></run></mzML>\n<indexList count="1"><index name="spectrum"><offset idRef="a">10</offset></index></indexList>\n'
            '<indexListOffset>123</indexListOffset>\n<fileChecksum>abc</fileChecksum>\n</indexedmzML>\n')


VELOS = mzml([
    spectrum(1, 1, 'FTMS + p ESI Full ms [200.00-2000.00]', None),
    spectrum(2, 2, 'FTMS + c ESI d Full ms2 664.46@hcd30.00 [90.00-675.00]', 'beam-type collision-induced dissociation'),
    spectrum(3, 2, 'ITMS + c ESI d Full ms2 664.46@cid35.00 [170.00-675.00]', 'collision-induced dissociation'),
    spectrum(4, 1, 'FTMS + p ESI Full ms [200.00-2000.00]', None),
])


def test_keeps_ms1_and_hcd_and_drops_cid():
    out, counts = keep_activation(VELOS, 'hcd')
    assert counts == {'ms1': 2, 'kept_msn': 1, 'dropped_msn': 1}
    assert 'scan=2' in out and 'scan=3' not in out
    assert out.count('<spectrum ') == 3


def test_output_is_plain_consistent_mzml():
    out, _ = keep_activation(VELOS, 'hcd')
    assert '<indexedmzML' not in out and '<indexList' not in out and '<fileChecksum>' not in out
    assert re.search(r'<spectrumList count="3"', out)
    assert [int(i) for i in re.findall(r'<spectrum [^>]*index="(\d+)"', out)] == [0, 1, 2]
    assert 'startTimeStamp="2020-02-19T03:14:17Z"' in out
    assert out.rstrip().endswith('</mzML>')


def test_uses_activation_term_when_filter_string_is_missing():
    data = mzml([
        spectrum(1, 1, None, None),
        spectrum(2, 2, None, 'beam-type collision-induced dissociation'),
        spectrum(3, 2, None, 'collision-induced dissociation'),
    ])
    out, counts = keep_activation(data, 'hcd')
    assert counts == {'ms1': 1, 'kept_msn': 1, 'dropped_msn': 1}
    assert 'scan=2' in out and 'scan=3' not in out
    out_cid, counts_cid = keep_activation(data, 'cid')
    assert counts_cid == {'ms1': 1, 'kept_msn': 1, 'dropped_msn': 1}
    assert 'scan=3' in out_cid and 'scan=2' not in out_cid


def test_file_without_the_activation_keeps_nothing_but_ms1():
    data = mzml([spectrum(1, 1, 'FTMS + p ESI Full ms', None), spectrum(2, 2, 'ITMS ms2 500.00@cid35.00', 'collision-induced dissociation')])
    out, counts = keep_activation(data, 'hcd')
    assert counts == {'ms1': 1, 'kept_msn': 0, 'dropped_msn': 1}
    assert out.count('<spectrum ') == 1
