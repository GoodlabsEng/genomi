"""Derived identity must retain the source binding and never overwrite original calls."""
import sqlite3
import unittest
from genomi.active_genome_index.array_genotypes import classify_array_genotype, called_genotype_tokens
from genomi.active_genome_index.normalized_alleles import write_normalized_allele_index, normalized_allele_index_metadata
from genomi.active_genome_index.allele_resolver import assess_call, assess_overlapping_calls

class NormalizedAlleleIndexTests(unittest.TestCase):

    def setUp(self):
        self.db = sqlite3.connect(':memory:')
        self.db.execute('create table metadata(key text primary key,value text)')
        self.db.execute('create table records(chrom text,pos integer,record_kind text,ref text,alt text,genotype text,observed_alleles text)')
        self.db.execute('insert into records values(\'1\',152,\'variant_call\',\'AA\',\'A\',\'1|1\',\'["A","A"]\')')
        self.original = self.db.execute('select * from records').fetchall()
        self.db.commit()
        self.metadata = dict(schemaVersion=1, genomeBuild='GRCh37', normalizationVersion='bcftools1.16', referenceSha256='a' * 64, sourceSha256='b' * 64, excludedContigs=[], observedAlleleCount=1)
        self.row = ('1', 'chr1', 150, 'GA', 'G', 1, 'A', 1, 'exact_allele', 'variant_call')

    def tearDown(self):
        self.db.close()

    def write(self, rows=None, metadata=None):
        return write_normalized_allele_index(self.db, expected_aliases={'1': 'chr1', 'chr1': 'chr1'}, normalized_rows=rows if rows is not None else [self.row], metadata=metadata if metadata is not None else self.metadata)

    def test_shifted_identity_preserves_original_row_and_is_complete(self):
        result = self.write()
        self.assertTrue(result['complete'])
        self.assertEqual(result['observedAlleleCount'], 1)
        self.assertEqual(normalized_allele_index_metadata(self.db), result)
        self.assertEqual(self.db.execute('select * from records').fetchall(), self.original)
        self.assertEqual(self.db.execute('select pos,ref,alt,record_rowid,original_allele_index from normalized_alleles').fetchone(), (150, 'GA', 'G', 1, 1))

    def test_invalid_binding_rolls_back_and_retry_can_complete(self):
        with self.assertRaisesRegex(ValueError, 'source_missing'):
            self.write([(*self.row[:5], 2, *self.row[6:])])
        self.assertIsNone(self.db.execute("select 1 from sqlite_schema where name='normalized_alleles'").fetchone())
        self.assertTrue(self.write()['complete'])

    def test_wrong_count_and_duplicate_initialization_are_rejected(self):
        with self.assertRaisesRegex(ValueError, 'count_mismatch'):
            self.write(metadata={**self.metadata, 'observedAlleleCount': 2})
        self.write()
        with self.assertRaisesRegex(ValueError, 'already_exists'):
            self.write(metadata={**self.metadata, 'sourceSha256': 'c' * 64})

    def test_zero_observed_alleles_is_complete_not_a_clinical_negative(self):
        result = self.write([], {**self.metadata, 'observedAlleleCount': 0})
        self.assertEqual(result['observedAlleleCount'], 0)
        self.assertEqual(self.db.execute('select * from records').fetchall(), self.original)

class SourceAssessmentTests(unittest.TestCase):

    def record(self, genotype='0|1', depth='30', quality='60', alt='C'):
        return dict(chrom='1', pos=10, end=10, ref='A', alt=alt, genotype=genotype, filter='PASS', is_variant=True, sample={'GT': genotype, 'DP': depth, 'GQ': quality, 'PS': '10'}, record_kind='variant_call')

    def test_phase_haploid_partial_and_unreported_alt(self):
        phase = assess_call(self.record(alt='C,G'))
        self.assertEqual(phase['alleles'], ['A', 'C'])
        self.assertTrue(phase['phased'])
        self.assertEqual(phase['phase_set'], '10')
        self.assertEqual(assess_call(self.record('1'))['ploidy'], 1)
        self.assertEqual(assess_call(self.record('1/.'))['status'], 'no_call')

    def test_duplicate_confidence_is_weakest_and_conflicts_are_not_selected(self):
        self.assertEqual(assess_overlapping_calls([self.record(), self.record(depth='4')])['quality'], 'low')
        self.assertEqual(assess_overlapping_calls([self.record(), self.record('1|1')])['status'], 'conflicting')

    def test_reordered_alleles_have_same_copies_and_phase_ambiguity_is_preserved(self):
        other = self.record('0|2', alt='G,C')
        self.assertEqual(assess_overlapping_calls([self.record(), other])['status'], 'called')
        reversed_phase = self.record('1|0')
        result = assess_overlapping_calls([self.record(), reversed_phase])
        self.assertEqual(result['source_call_count'], 2)
        self.assertIsNone(result['phased'])
        self.assertEqual([c['raw_genotype'] for c in result['source_calls']], ['0|1', '1|0'])
class ArrayGenotypeClassificationTests(unittest.TestCase):
    def test_supported_pairs_preserve_original_and_tokens(self):
        for genotype, alleles in (("A/G", ("A", "G")), ("C/C", ("C", "C")),
                                  ("AG", ("A", "G")), ("c/g", ("C", "G"))):
            with self.subTest(genotype=genotype):
                call = classify_array_genotype(genotype)
                self.assertEqual(call.raw_genotype, genotype)
                self.assertEqual(call.alleles, alleles)
                self.assertEqual(call.status, "called")
                self.assertTrue(call.snp_eligible)

    def test_missing_tokens_are_not_complete_haploid_calls(self):
        for genotype in ("./A", "A/.", "./.", "N/N", "00/A", "NN/A", "A/-", "--/--", "A0", "AN", "A-", "--", ""):
            with self.subTest(genotype=genotype):
                call = classify_array_genotype(genotype)
                self.assertEqual(call.status, "no_call")
                self.assertFalse(call.snp_eligible)
                self.assertEqual(called_genotype_tokens(genotype), [])
                assessed = assess_call({"record_kind": "array_no_call", "genotype": genotype}, array=True)
                self.assertEqual(assessed["status"], "no_call")
                self.assertEqual(assessed["alleles"], [])
        self.assertEqual(classify_array_genotype("./A").alleles, (".", "A"))

    def test_complex_symbolic_and_legacy_indel_tokens_cannot_establish_snp_evidence(self):
        for genotype in ("A/AT", "AT/A", "AT/AT", "B/B", "XYZ/A", "A/<DEL>", "*/A", "I/D", "ID", "A/G/T"):
            with self.subTest(genotype=genotype):
                call = classify_array_genotype(genotype)
                self.assertFalse(call.snp_eligible)
                assessed = assess_call({"record_kind": "array_call", "genotype": genotype}, array=True)
                self.assertEqual(assessed["status"], "unresolved")
                self.assertEqual(assessed["raw_genotype"], genotype)
        self.assertEqual(classify_array_genotype("A/AT").alleles, ("A", "AT"))
        self.assertEqual(classify_array_genotype("ID").status, "called")

    def test_numeric_vcf_helpers_preserve_complete_calls_and_reject_partial(self):
        for genotype, tokens in (("0|1", ["0", "1"]), ("1", ["1"]), ("0/2", ["0", "2"]),
                                 ("1/.", []), ("./0", [])):
            self.assertEqual(called_genotype_tokens(genotype), tokens)


    def test_normalized_writer_rejects_forged_partial_or_complex_components(self):
        for genotype, observed in (("A/AT", '["A", "AT"]'), ("AT/A", '["AT", "A"]'),
                                   ("./A", '["A"]'), ("A/G", '["A", "AT"]')):
            with self.subTest(genotype=genotype), sqlite3.connect(":memory:") as db:
                db.executescript("create table metadata(key text primary key,value text);"
                                 "create table records(chrom text,pos integer,record_kind text,ref text,alt text,genotype text,observed_alleles text);")
                db.execute("insert into records values('1',10,'array_call','.','.',?,?)", (genotype, observed))
                with self.assertRaisesRegex(ValueError, "source_allele_mismatch"):
                    write_normalized_allele_index(db, expected_aliases={"1": "chr1"},
                                                  normalized_rows=[("1:A:array", "chr1", 10, ".", "A", 1, "A", None, "array_allele", "array_call")],
                                                  metadata={"schemaVersion": 1, "genomeBuild": "GRCh37", "normalizationVersion": "synthetic",
                                                            "referenceSha256": "a" * 64, "sourceSha256": "b" * 64,
                                                            "excludedContigs": [], "observedAlleleCount": 1})
                self.assertIsNone(db.execute("select 1 from metadata where key='normalized_allele_index'").fetchone())

if __name__ == '__main__':
    unittest.main()
