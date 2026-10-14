"""Synthetic BGZF, original-call and batch-resolver regression coverage."""
from pathlib import Path
import hashlib, json, sqlite3, tempfile, time, unittest
from unittest.mock import patch
import pysam
from pysam.libcbgzf import BGZFile
from genomi.active_genome_index import _agi_schema as schema
from genomi.active_genome_index.reader import ActiveGenomeIndexReader, ActiveGenomeIndexNeed
from genomi.active_genome_index.vcf import parse_record_fields
from genomi.active_genome_index.source_intake.agi_store import _array_record_row
from genomi.active_genome_index.normalized_alleles import write_normalized_allele_index
from genomi.active_genome_index import allele_resolver as resolver
OUT = None
CASES = {}

def source_row(pos, ref='A', alt='C', gt='0/1', fmt='GT:DP:GQ', sample=None, info='.', identifier='.'):
    if sample is None:
        sample = gt + ':30:99'
    return f'1\t{pos}\t{identifier}\t{ref}\t{alt}\t60\tPASS\t{info}\t{fmt}\t{sample}\n'

def fixture(name, lines, keys=(), provider='vcf', array=None, reference_overrides=()):
    root = OUT / name
    root.mkdir(exist_ok=True)
    reference = root / 'reference.fa'
    bases = list('A' * 1000)
    bases[99:102] = list('GAT')
    for pos, base in reference_overrides:
        bases[pos - 1] = base
    reference.write_text('>chr1\n' + ''.join(bases) + '\n')
    pysam.faidx(str(reference))
    digest = hashlib.sha256(reference.read_bytes()).hexdigest()
    index = root / 'active-genome-index.sqlite'
    if index.exists():
        index.unlink()
    db = sqlite3.connect(index)
    schema._reset_schema(db)
    source = root / 'source.vcf.gz'
    records = []
    if array is None:
        with BGZFile(str(source), 'wb') as handle:
            handle.write(b'##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTHETIC\n')
            for line in lines:
                raw = line.encode()
                offset = handle.tell()
                handle.write(raw)
                record = parse_record_fields(line.rstrip().split('\t'), sample_name='SYNTHETIC', sample_index=0, offset=offset, line_length=len(raw))
                records.append(schema._record_row(record))
        Path(str(source) + '.gzi').touch()
    else:
        for (i, (identifier, pos, gt)) in enumerate(array, 1):
            called = bool(gt) and gt not in ('--', '00', 'NN', '.')
            records.append(_array_record_row(dict(rsid=identifier, chrom='1', pos=pos, genotype=gt), row_index=i, is_called=called, sample_name='SYNTHETIC', source_format=provider))
    schema._insert_record_batch(db, records)
    schema._create_query_indexes(db)
    metadata = dict(schema_version=schema.SCHEMA_VERSION, source_format=provider, genome_build='GRCh37', active_genome_index_complete=True, active_genome_index_build_status='completed')
    if array is None:
        stat = source.stat()
        metadata.update(vcf_path=str(source), vcf_size_bytes=stat.st_size, vcf_mtime_ns=stat.st_mtime_ns)
    db.executemany('insert into metadata values (?,?)', [(key, json.dumps(value)) for (key, value) in metadata.items()])
    normalized = []
    for (rowid, pos, ref, alt, original, allele_index) in keys:
        kind = 'array_call' if array is not None else 'variant_call'
        normalized.append((f'{rowid}:{original}', 'chr1', pos, ref, alt, rowid, original, allele_index, 'exact_allele', kind))
    if array is not None:
        for (rowid, (_, pos, gt)) in enumerate(array, 1):
            for base in sorted(set(gt) & set('ACGT')):
                normalized.append((f'{rowid}:{base}:array', 'chr1', pos, '.', base, rowid, base, None, 'array_allele', 'array_call'))
    write_normalized_allele_index(db, expected_aliases={'1': 'chr1', 'chr1': 'chr1'}, normalized_rows=normalized, metadata=dict(schemaVersion=1, genomeBuild='GRCh37', normalizationVersion='synthetic_known_keys', referenceSha256=digest, sourceSha256='b' * 64, excludedContigs=[], observedAlleleCount=len(normalized)))
    db.commit()
    original = db.execute('select * from records').fetchall()
    db.close()
    reader = ActiveGenomeIndexReader(index, ActiveGenomeIndexNeed.NONE, readiness={'complete': True}, genome_build='GRCh37')

    def resolve(targets, observed_only=False, deadline=None):
        result = reader.resolve_allele_targets(targets, reference_fasta=reference, expected_reference_sha256=digest, deadline=deadline or time.monotonic() + 10, observed_only=observed_only)
        with sqlite3.connect(index) as current:
            assert current.execute('select * from records').fetchall() == original, 'original records changed'
        return result
    return (resolve, reader, reference, digest)

def target(pos, ref='A', alt='C', rsid=None):
    result = dict(chrom='1', pos=pos, ref=ref, alt=alt)
    if rsid:
        result['rsid'] = rsid
    return result

def run(name, check):
    CASES[name] = check

def status(row, status, exact=None, dosage='skip', assessment=None):
    assert row['status'] == status, row
    if exact is not None:
        assert row['exact_match'] is exact, row
    if dosage != 'skip':
        assert row['target_dosage'] == dosage, row
    if assessment:
        assert row['call_assessment']['status'] == assessment, row
    return {key: row.get(key) for key in ('status', 'exact_match', 'target_dosage', 'reference_reason', 'reference_scope')}

def padded():
    (resolve, *_) = fixture('padded', [source_row(100, 'GAT', 'GCT', gt='0|1', fmt='GT:DP:GQ:PS', sample='0|1:30:99:42')], [(1, 101, 'A', 'C', 'GCT', 1)])
    row = resolve([target(101)])['results'][0]
    assert row['call_assessment']['phased'] and row['call_assessment']['phase_set'] == '42', row
    assert row['source_records'][0]['pos'] == 100 and row['source_records'][0]['ref'] == 'GAT', row
    return status(row, 'observed', True, 1, 'called')
run('padded_exact_and_source_phase', padded)

def projected():
    (resolve, *_) = fixture('projected', [source_row(100, 'GAT', 'GCA')], [(1, 101, 'AT', 'CA', 'GCA', 1)])
    row = resolve([target(101)])['results'][0]
    assert row['site_projection']['allele_bases'] == ['A', 'C'], row
    return status(row, 'unavailable', False, None, 'called')
run('larger_GCA_is_projection_only', projected)

def shifted():
    (resolve, *_) = fixture('shifted', [source_row(150, 'AA', 'A'), source_row(151, 'A', 'C')], [(1, 103, 'AA', 'A', 'A', 1), (2, 151, 'A', 'C', 'C', 1)])
    row = resolve([target(103, 'AA', 'A')])['results'][0]
    assert {record['pos'] for record in row['source_records']} == {150, 151}, row
    return status(row, 'unavailable', True, None, 'unresolved')
run('shifted_deletion_interior_overlap', shifted)

def multiallelic():
    (resolve, *_) = fixture('multiallelic', [source_row(200, alt='C,G', gt='0/2'), source_row(201, gt='1'), source_row(202, gt='1/.')], [(1, 200, 'A', 'G', 'G', 2), (2, 201, 'A', 'C', 'C', 1)])
    rows = resolve([target(200), target(200, alt='G'), target(201), target(202)])['results']
    status(rows[0], 'not_observed', False, 0, 'called')
    status(rows[1], 'observed', True, 1, 'called')
    status(rows[2], 'observed', True, 1, 'called')
    assert rows[2]['call_assessment']['ploidy'] == 1, rows[2]
    status(rows[3], 'no_call', False, None, 'no_call')
    return [row['status'] for row in rows]
run('unreported_ALT_haploid_partial', multiallelic)

def phase_duplicate():
    (resolve, *_) = fixture('phase', [source_row(210, gt='0|1', fmt='GT:DP:GQ:PS', sample='0|1:30:99:42'), source_row(210, gt='1|0', fmt='GT:DP:GQ:PS', sample='1|0:4:99:42')], [(1, 210, 'A', 'C', 'C', 1), (2, 210, 'A', 'C', 'C', 1)])
    row = resolve([target(210)])['results'][0]
    assert row['call_assessment']['phased'] is None and row['call_assessment']['phase_set'] is None, row
    assert row['call_assessment']['quality'] == 'low' and row['call_assessment']['depth'] == 4, row
    return status(row, 'observed', True, 1, 'called')
run('duplicate_phase_and_weakest_quality', phase_duplicate)

def vcf_absent():
    (resolve, *_) = fixture('absent', [source_row(200)], [(1, 200, 'A', 'C', 'C', 1)])
    return status(resolve([target(400)])['results'][0], 'not_observed', False, None, 'unresolved')
run('plain_VCF_absence_inconclusive', vcf_absent)

def block(name, sample='0/0:30:99:30', alt='<NON_REF>', fmt='GT:DP:GQ:MIN_DP', extra=()):
    return fixture(name, [source_row(300, alt=alt, fmt=fmt, sample=sample, info='END=310'), *extra], provider='gvcf')[0]
run('gVCF_valid_reference', lambda : status(block('block')([target(305)])['results'][0], 'reference_confirmed', False, 0, 'called'))
run('gVCF_low_MIN_DP', lambda : status(block('block_low', sample='0/0:30:99:3')([target(305)])['results'][0], 'unavailable', False, None, 'called'))
run('gVCF_missing_MIN_DP', lambda : status(block('block_missing', sample='0/0:30:99', fmt='GT:DP:GQ')([target(305)])['results'][0], 'unavailable', False, None, 'called'))
run('gVCF_invalid_LEN_span', lambda : status(block('block_len', sample='0/0:30:99:30:2', fmt='GT:DP:GQ:MIN_DP:LEN')([target(305)])['results'][0], 'unavailable', False, None, 'called'))
run('gVCF_structural_END_is_not_confidence', lambda : status(block('block_structural', alt='<DEL>')([target(305)])['results'][0], 'unavailable', False, None, 'called'))
run('gVCF_overlapping_reference_rows', lambda : status(block('block_overlap', extra=[source_row(305, alt='<NON_REF>', fmt='GT:DP:GQ:MIN_DP', sample='0/0:30:99:30', info='END=310')])([target(305)])['results'][0], 'unavailable', False, None, 'unresolved'))

def arr(name, array, expected, assessment, exact=False, dosage=None):
    (resolve, *_) = fixture(name, [], provider='23andme', array=array)
    row = resolve([target(101, rsid='rs1004')])['results'][0]
    if assessment == 'called':
        assert row['call_assessment']['quality'] == 'unknown', row
    if expected == 'reference_confirmed':
        assert row['reference_scope'] == 'assayed_array_marker', row
    return status(row, expected, exact, dosage, assessment)
run('array_reference', lambda : arr('array_ref', [('rs1004', 101, 'AA')], 'reference_confirmed', 'called', False, 0))
run('array_missing', lambda : arr('array_missing', [('rs1004', 101, '--')], 'no_call', 'no_call'))
run('array_unrelated_missing_probe', lambda : arr('array_other', [('rs1004', 101, 'AC'), ('i100', 101, '--')], 'observed', 'called', True, 1))
run('array_same_probe_missing', lambda : arr('array_same', [('rs1004', 101, 'AC'), ('rs1004', 101, '--')], 'no_call', 'no_call', True, None))
run('array_conflicting_probes', lambda : arr('array_conflict', [('rs1004', 101, 'AC'), ('i100', 101, 'CC')], 'unavailable', 'conflicting', True, None))

def array_unknown_target(name, ref, alt):
    (resolve, *_) = fixture(name, [], provider='23andme', array=[('rs1004', 101, 'CC')],
                          reference_overrides=[(101, ref)])
    envelope = resolve([target(101, ref=ref, alt=alt, rsid='rs1004')])
    row = envelope['results'][0]
    assert envelope['relevant_target_ordinals'] == [], envelope
    assert row['reference_reason'] == 'array_allele_model_unsupported', row
    return status(row, 'unavailable', False, None, 'unresolved')
run('array_unknown_reference_does_not_invent_absence', lambda : array_unknown_target('array_unknown_ref', 'N', 'C'))
run('array_unknown_alternate_is_unsupported', lambda : array_unknown_target('array_unknown_alt', 'A', 'N'))

def array_negative_requires_target_alleles():
    for genotype in ('TC', 'AC'):
        name = f'array_negative_{genotype}'
        (resolve, *_) = fixture(name, [], provider='23andme',
                               array=[('rs1004', 200, genotype)])
        row = resolve([target(200, alt='G', rsid='rs1004')])['results'][0]
        assert row['reference_reason'] == 'array_allele_model_unsupported', row
        assert row['call_assessment']['alleles'] == list(genotype), row
        assert row['call_assessment']['quality'] == 'unknown', row
        assert row['source_records'][0]['genotype'] == genotype, row
        assert 'reference_scope' not in row, row
        status(row, 'unavailable', False, None, 'called')
run('array_negative_requires_target_alleles', array_negative_requires_target_alleles)

def array_literal_positives_keep_third_alleles():
    for genotype in ('AG', 'CG', 'TG'):
        name = f'array_positive_{genotype}'
        (resolve, *_) = fixture(name, [], provider='23andme',
                               array=[('rs1004', 200, genotype)])
        row = resolve([target(200, alt='G', rsid='rs1004')])['results'][0]
        assert row['call_assessment']['alleles'] == list(genotype), row
        assert row['call_assessment']['quality'] == 'unknown', row
        assert row['source_records'][0]['genotype'] == genotype, row
        status(row, 'observed', True, 1, 'called')
run('array_literal_positives_keep_third_alleles', array_literal_positives_keep_third_alleles)

def other_snv_cannot_supply_indel_zero():
    for provider in ('vcf', 'gvcf'):
        (resolve, *_) = fixture(f'other_indel_{provider}', [source_row(200)],
                               [(1, 200, 'A', 'C', 'C', 1)], provider=provider)
        row = resolve([target(200, alt='AT')])['results'][0]
        assert row['reference_reason'] == 'sequence_allele_not_identical', row
        assert row['call_assessment']['alleles'] == ['A', 'C'], row
        assert row['call_assessment']['quality'] == 'pass', row
        assert row['source_records'][0]['genotype'] == '0/1', row
        status(row, 'unavailable', False, None, 'called')
run('other_snv_cannot_supply_indel_zero', other_snv_cannot_supply_indel_zero)

def other_snv_zero_requires_qualified_call():
    cases = (
        ('pass', 'PASS', 'GT:DP:GQ', '0/1:30:99', 'pass', None),
        ('low', 'PASS', 'GT:DP:GQ', '0/1:4:5', 'low', 'depth_below_threshold'),
        ('unknown', 'PASS', 'GT:DP:GQ', '0/1:.:.', 'unknown', 'depth_missing'),
        ('sample_filter', 'PASS', 'GT:DP:GQ:FT', '0/1:30:99:LowQual', 'low', 'sample_filter_failed'),
        ('record_filter', 'LowQual', 'GT:DP:GQ', '0/1:30:99', 'low', 'record_filter_failed'),
        ('unapplied_filter', '.', 'GT:DP:GQ', '0/1:30:99', 'unknown', 'record_filter_not_applied'),
    )
    for provider in ('vcf', 'gvcf'):
        for name, record_filter, fmt, sample, quality, reason in cases:
            line = source_row(200, fmt=fmt, sample=sample).replace('\tPASS\t', f'\t{record_filter}\t')
            (resolve, *_) = fixture(f'other_snv_{provider}_{name}', [line],
                                   [(1, 200, 'A', 'C', 'C', 1)], provider=provider)
            absent, observed = resolve([target(200, alt='G'), target(200)])['results']
            assert absent['call_assessment']['alleles'] == ['A', 'C'], absent
            assert absent['call_assessment']['quality'] == quality, absent
            assert absent['source_records'][0]['genotype'] == '0/1', absent
            assert absent['reference_reason'] == reason, absent
            if reason is None:
                status(absent, 'not_observed', False, 0, 'called')
            else:
                assert reason in absent['call_assessment']['quality_reasons'], absent
                status(absent, 'unavailable', False, None, 'called')
            assert observed['call_assessment']['quality'] == quality, observed
            status(observed, 'observed', True, 1, 'called')
run('other_snv_zero_requires_qualified_call', other_snv_zero_requires_qualified_call)

def multi_ref_qualified_anchor():
    for provider in ('vcf', 'gvcf'):
        for alt in ('G', 'G,<NON_REF>'):
            (resolve, *_) = fixture(f'anchor_{provider}_{alt}', [source_row(
                100, ref='GAT', alt=alt, fmt='GT:DP:GQ:MIN_DP', sample='0/0:30:99:30',
            )], provider=provider)
            row = resolve([target(100, ref='G', alt='C')])['results'][0]
            assert row['call_assessment']['alleles'] == ['GAT', 'GAT'], row
            assert row['source_records'][0]['ref'] == 'GAT', row
            assert row['source_records'][0]['genotype'] == '0/0', row
            status(row, 'reference_confirmed', False, 0, 'called')
            interior = resolve([target(101)])['results'][0]
            reason = ('plain_vcf_reference_projection_unsupported' if provider == 'vcf'
                      else 'unsupported_source_reference_span')
            assert interior['reference_reason'] == reason, interior
            status(interior, 'unavailable', False, None, 'called')
            indel = resolve([target(100, ref='G', alt='GT')])['results'][0]
            assert indel['reference_reason'] == 'unsupported_reference_projection', indel
            status(indel, 'unavailable', False, None, 'called')
run('multi_ref_qualified_anchor_and_projection_bounds', multi_ref_qualified_anchor)

def multi_ref_anchor_vetoes():
    cases = (
        ('wrong_ref', 'GAC', '0/0:30:99:30', (), 'unavailable', 'source_reference_sequence_mismatch'),
        ('overlap', 'GAT', '0/0:30:99:30', (source_row(100, ref='G', alt='C'),),
         'unavailable', 'source_genotype_unresolved'),
        ('partial', 'GAT', '0/.:30:99:30', (), 'no_call', 'source_genotype_no_call'),
        ('low', 'GAT', '0/0:4:5:4', (), 'unavailable', 'depth_below_threshold'),
        ('unknown', 'GAT', '0/0:30:99:.', (), 'unavailable', 'minimum_depth_missing'),
    )
    for provider in ('vcf', 'gvcf'):
        for name, ref, sample, extra, expected, reason in cases:
            (resolve, *_) = fixture(f'anchor_veto_{provider}_{name}', [source_row(
                100, ref=ref, alt='G,<NON_REF>', fmt='GT:DP:GQ:MIN_DP', sample=sample,
            ), *extra], provider=provider)
            row = resolve([target(100, ref='G', alt='C')])['results'][0]
            assert row['reference_reason'] == reason, row
            status(row, expected, False, None)
run('multi_ref_anchor_vetoes', multi_ref_anchor_vetoes)

def timeout_partial():
    (resolve, *_) = fixture('timeout', [source_row(200), source_row(210)], [(1, 200, 'A', 'C', 'C', 1), (2, 210, 'A', 'C', 'C', 1)])
    original = resolver._resolve_target
    monotonic = time.monotonic
    deadline = monotonic() + 10
    expired = [False]
    calls = [0]

    def clock():
        return deadline + 1 if expired[0] else monotonic()

    def interrupted(*args, **kwargs):
        calls[0] += 1
        if calls[0] == 2:
            expired[0] = True
            raise sqlite3.OperationalError('interrupted')
        return original(*args, **kwargs)
    with patch.object(resolver, '_resolve_target', interrupted), patch.object(resolver.time, 'monotonic', clock):
        envelope = resolve([target(200), target(210)], deadline=deadline)
    assert envelope['status'] == 'partial' and envelope['returned_target_count'] == 1 and (envelope['omitted_target_count'] == 1), envelope
    assert envelope['relevant_target_count'] == 2 and envelope['relevant_target_ordinals'] == [0, 1], envelope
    return {key: envelope[key] for key in ('status', 'returned_target_count', 'omitted_target_count', 'unavailable_reason')}
run('SQLite_timeout_preserves_prior_rows', timeout_partial)

def symbolic_unreported_span():
    (resolve, *_) = fixture('symbolic_unreported', [source_row(400, alt='C,<DEL>', info='END=450'), source_row(425, alt='G')], [(1, 400, 'A', 'C', 'C', 1), (2, 425, 'A', 'G', 'G', 1)])
    row = resolve([target(400)])['results'][0]
    return status(row, 'observed', True, 1, 'called')
run('unreported_symbolic_END_does_not_extend_called_SNP', symbolic_unreported_span)

def inventory_cap():
    lines = [source_row(pos) for pos in range(200, 301)]
    keys = [(i, pos, 'A', 'C', 'C', 1) for (i, pos) in enumerate(range(200, 301), 1)]
    (resolve, *_) = fixture('inventory_cap', lines, keys)
    envelope = resolve([target(pos) for pos in range(200, 10200)], observed_only=True)
    assert envelope['total_target_count'] == 10000 and envelope['relevant_target_count'] == 101, envelope
    assert envelope['returned_target_count'] == 100 and envelope['omitted_target_count'] == 1, envelope
    assert envelope['relevant_target_ordinals'] == list(range(101)), envelope
    assert envelope['status'] == 'partial' and envelope['unavailable_reason'] == 'resolution_target_limit', envelope
    assert all((row['status'] == 'observed' for row in envelope['results'])), envelope
    return {key: envelope[key] for key in ('total_target_count', 'relevant_target_count', 'returned_target_count', 'omitted_target_count', 'status')}
run('10k_candidate_inventory_assessment_cap', inventory_cap)

def source_overflow():
    lines = [source_row(200) for _ in range(65)] + [source_row(210)]
    keys = [(i, 200, 'A', 'C', 'C', 1) for i in range(1, 66)] + [(66, 210, 'A', 'C', 'C', 1)]
    (resolve, *_) = fixture('source_overflow', lines, keys)
    envelope = resolve([target(200), target(210)], observed_only=True)
    assert envelope['status'] == 'complete' and envelope['relevant_target_count'] == 2, envelope
    status(envelope['results'][0], 'unavailable', True, None, 'unresolved')
    assert envelope['results'][0]['reference_reason'] == 'source_call_limit', envelope
    status(envelope['results'][1], 'observed', True, 1, 'called')
    return [row['status'] for row in envelope['results']]
run('source_overflow_is_target_local', source_overflow)

def no_inferred_confidence():
    line = source_row(200, fmt='GT:AD:PL', sample='0/1:15,15:99,0,99', info='DP=30')
    (resolve, *_) = fixture('no_inferred_metrics', [line], [(1, 200, 'A', 'C', 'C', 1)])
    row = resolve([target(200)])['results'][0]
    assert row['call_assessment']['quality'] == 'unknown', row
    assert row['call_assessment']['depth'] is None and row['call_assessment']['genotype_quality'] is None, row
    assert row['source_records'][0]['depth'] is None and row['source_records'][0]['genotype_quality'] is None, row
    return status(row, 'observed', True, 1, 'called')
run('INFO_AD_PL_do_not_infer_sequencing_confidence', no_inferred_confidence)

def empty_metadata_gate():
    (resolve, reader, *_) = fixture('empty_gate', [])
    envelope = resolve([])
    assert envelope['status'] == 'complete' and envelope['total_target_count'] == 0, envelope
    with sqlite3.connect(reader.agi_path) as db:
        raw = json.loads(db.execute("select value from metadata where key='normalized_allele_index'").fetchone()[0])
        raw['schemaVersion'] = 0
        db.execute("update metadata set value=? where key='normalized_allele_index'", (json.dumps(raw),))
    try:
        resolve([])
    except ValueError as exc:
        assert str(exc) == 'normalized_allele_index_incomplete', str(exc)
        return 'empty invalid index refused'
    raise AssertionError('empty invalid index accepted')
run('empty_request_requires_compatible_enrichment', empty_metadata_gate)

class NativeAlleleResolverTests(unittest.TestCase):

    def setUp(self):
        global OUT
        self.temporary = tempfile.TemporaryDirectory(prefix='genomi-resolver-test-')
        OUT = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

def _test_case(check):

    def test(self):
        check()
    return test
for (_name, _check) in CASES.items():
    setattr(NativeAlleleResolverTests, 'test_' + _name, _test_case(_check))
if __name__ == '__main__':
    unittest.main()
