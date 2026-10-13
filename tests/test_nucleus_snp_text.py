"""Synthetic headerless SNP text through validation, original observations and exact resolution."""
from __future__ import annotations

import hashlib
import json
import os
import sqlite3
import tempfile
import time
import unittest
from pathlib import Path
from unittest.mock import patch

import pysam

from genomi.active_genome_index.array_genotypes import classify_array_genotype, called_genotype_tokens
from genomi.active_genome_index.allele_resolver import assess_call
from genomi.active_genome_index.normalized_alleles import write_normalized_allele_index
from genomi.active_genome_index.reader import ActiveGenomeIndexReader, ActiveGenomeIndexNeed
from genomi.active_genome_index.source_intake.detection import detect_source
from genomi.active_genome_index.source_intake.dispatch import parse_source
from genomi.active_genome_index.source_intake.nucleus_text import MARKER_DEFINITIONS


def compatible_rows() -> list[str]:
    # Public coordinates with entirely synthetic genotypes, not a personal sample.
    genotypes = ("A/A", "G/G", "G/G", "T/T", "T/T", "C/T", "G/G", "A/T")
    return [f'{marker["rsid"]}\t{marker["chromosome"]}\t{marker["grch37Position"]}\t{genotype}\n'
            for marker, genotype in zip(MARKER_DEFINITIONS["markers"], genotypes)]


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


class NucleusSnpTextTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.environment = patch.dict(os.environ, {"GENOMI_HOME": str(self.root / "genomi"),
                                                  "GENOMI_SESSION_ID": "", "GENOMI_CONTEXT": ""})
        self.environment.start()
        self.addCleanup(self.environment.stop)
        self.source = self.root / "synthetic.txt"

    def write(self, rows=None):
        self.source.write_text("".join(compatible_rows() if rows is None else rows), encoding="utf-8")
        return self.source

    def test_headerless_detection_is_validated_not_vendor_authenticated(self):
        detection = detect_source(self.write())
        self.assertEqual(detection.source_format, "nucleus")
        self.assertEqual(detection.array_dialect, "nucleus_snp_text")
        self.assertEqual(detection.reference_build, "GRCh37")
        self.assertIsNone(detection.provider)
        self.assertEqual(detection.reference_validation,
                         {"markerSetVersion": 1, "coordinateMarkerCount": 8, "orientationMarkerCount": 5})
        self.assertEqual(detection.array_call_summary["total_rows"], 8)

    def test_cached_array_requires_complete_call_summary_and_validation_binding(self):
        self.write()
        first = parse_source(self.source)
        index = first["outputs"]["agi_path"]
        self.assertEqual(parse_source(self.source)["steps"][0]["result"]["status"], "cached")
        with sqlite3.connect(index) as db:
            db.execute("delete from metadata where key='array_call_summary'")
            db.commit()
        self.assertEqual(parse_source(self.source)["steps"][0]["result"]["status"], "completed")
        with sqlite3.connect(index) as db:
            db.execute("update metadata set value='null' where key='reference_validation'")
            db.commit()
        self.assertEqual(parse_source(self.source)["steps"][0]["result"]["status"], "completed")

    def test_alphabetic_unsupported_calls_are_retained_without_snp_keys(self):
        additions = ["rs900000008\tCHR1\t80\tB/B\n", "rs900000009\tchr1\t90\tXYZ/A\n"]
        self.write(compatible_rows() + additions)
        parsed = parse_source(self.source)
        index = Path(parsed["outputs"]["agi_path"])
        reader = ActiveGenomeIndexReader(index, ActiveGenomeIndexNeed.NONE, readiness={"complete": True}, genome_build="GRCh37")
        with sqlite3.connect(index) as db:
            db.row_factory = sqlite3.Row
            records = list(db.execute("select chrom,genotype,observed_alleles from records where pos in (80,90) order by pos"))
            self.assertEqual([tuple(row) for row in records], [("1", "B/B", None), ("1", "XYZ/A", None)])
            reader.stage_clinvar_match_records(db)
            self.assertEqual(db.execute("select count(*) from temp.selected_active_genome_index_records where pos in (80,90)").fetchone()[0], 0)
            summary = json.loads(db.execute("select value from metadata where key='array_call_summary'").fetchone()[0])
            self.assertEqual(summary["unsupported_calls"], 2)

    def test_full_file_validation_reads_markers_beyond_detection_probe(self):
        self.write(["rs900000001\t1\t10\tA/G\n"] * 15000 + compatible_rows())
        self.assertGreater(self.source.stat().st_size, 256 * 1024)
        self.assertEqual(detect_source(self.source).reference_validation["coordinateMarkerCount"], 8)

    def test_coordinate_and_orientation_contradictions_veto_other_valid_markers(self):
        additions = (("rs429358\t19\t44908684\tT/T\n", "coordinate_conflict"),
                     ("rs429358\t19\t45411941\tA/A\n", "orientation_conflict"),
                     ("rs762551\t15\t75041917\tA/T\n", "marker_alleles_unknown"))
        for extra, reason in additions:
            with self.subTest(reason=reason), self.assertRaisesRegex(ValueError, reason):
                detect_source(self.write(compatible_rows() + [extra]))

    def test_unknown_supported_base_calls_and_insufficient_evidence_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "evidence_insufficient"):
            detect_source(self.write(compatible_rows()[:3]))
        # All these synthetic calls are compatible with either orientation.
        rows = []
        for marker in MARKER_DEFINITIONS["markers"]:
            forward = set(marker["forwardAlleles"])
            reverse = {allele.translate(str.maketrans("ACGT", "TGCA")) for allele in forward}
            base = sorted(forward & reverse)[0] if forward & reverse else None
            genotype = f"{base}/{base}" if base else "./."
            rows.append(f'{marker["rsid"]}\t{marker["chromosome"]}\t{marker["grch37Position"]}\t{genotype}\n')
        with self.assertRaisesRegex(ValueError, "evidence_insufficient"):
            detect_source(self.write(rows))

    def test_duplicate_markers_and_single_chromosome_orientation_are_not_enough(self):
        with self.assertRaisesRegex(ValueError, "evidence_insufficient"):
            detect_source(self.write(compatible_rows()[:1] * 8))
        rows = compatible_rows()
        rows[0] = rows[0].replace("A/A", "C/C")
        rows[2] = rows[2].replace("G/G", "A/A")
        rows[3] = rows[3].replace("T/T", "C/C")
        with self.assertRaisesRegex(ValueError, "evidence_insufficient"):
            detect_source(self.write(rows))

    def test_malformed_rows_cannot_be_skipped(self):
        for row in ("rs900000001\t1\t10\tA/G\textra\n", "rs900000001\t1\t0\tA/G\n",
                    "rs900000001\t1\t10\tAG\n", "rs900000001\t25\t10\tA/G\n",
                    "# reference build 37\n", "rs900000001\t1\t10\tA//G\n"):
            with self.subTest(row_shape=row.split("\t")[0]), self.assertRaisesRegex(ValueError, "layout_invalid"):
                detect_source(self.write(compatible_rows() + [row]))

    def test_original_bytes_raw_genotype_summary_and_explicit_build_are_consistent(self):
        self.write(compatible_rows() + ["rs900000001\tchr1\t10\tc/g\n", "rs900000002\tM\t10\tA/A\n"])
        original = self.source.read_bytes()
        with self.assertRaisesRegex(ValueError, "build_conflict"):
            parse_source(self.source, genome_build="GRCh38")
        self.assertFalse((self.root / "genomi").exists())
        parsed = parse_source(self.source, genome_build="GRCh37")
        with sqlite3.connect(parsed["outputs"]["agi_path"]) as db:
            row = db.execute("select genotype,observed_alleles from records where rsid='rs900000001'").fetchone()
            metadata = {key: json.loads(value) for key, value in db.execute("select key,value from metadata")}
        self.assertEqual(row, ("c/g", '["C", "G"]'))
        self.assertEqual(metadata["array_dialect"], "nucleus_snp_text")
        self.assertEqual(metadata["reference_validation"]["markerSetVersion"], 1)
        self.assertEqual(metadata["array_call_summary"], {"total_rows": 10, "supported_single_base_calls": 10,
                                                       "incomplete_calls": 0, "unsupported_calls": 0})
        self.assertEqual(self.source.read_bytes(), original)
        self.assertIsNone(parsed["provider"])

    def _prepared_index(self):
        additions = ["rs900000001\t1\t10\tA/G\n", "rs900000002\t1\t20\tC/C\n",
                     "rs900000003\t1\t30\t./A\n", "rs900000004\t1\t40\tA/AT\n",
                     "rs900000005\t1\t50\tAT/A\n", "rs900000006\t1\t60\tA/<DEL>\n",
                     "rs900000007\t1\t70\tA/A\n"]
        self.write(compatible_rows() + additions)
        parsed = parse_source(self.source, genome_build="GRCh37")
        index = Path(parsed["outputs"]["agi_path"])
        reference = self.root / "reference.fa"
        reference.write_text(">chr1\n" + "A" * 1000 + "\n")
        pysam.faidx(str(reference))
        digest = hashlib.sha256(reference.read_bytes()).hexdigest()
        reader = ActiveGenomeIndexReader(index, ActiveGenomeIndexNeed.NONE, readiness={"complete": True}, genome_build="GRCh37")
        with sqlite3.connect(index) as db:
            db.row_factory = sqlite3.Row
            reader.stage_clinvar_match_records(db)
            aliases = {chrom: "chr" + chrom for chrom, in db.execute("select distinct chrom from records")}
            aliases.update({canonical: canonical for canonical in list(aliases.values())})
            selected = list(db.execute("select clinvar_batch_id,chrom,pos,ref,clinvar_match_alt,record_rowid,record_kind from temp.selected_active_genome_index_records"))
            rows = [(batch, aliases[chrom], pos, ref, alt, rowid, alt, None, "consumer_array_allele_inference", kind)
                    for batch, chrom, pos, ref, alt, rowid, kind in selected]
            write_normalized_allele_index(db, expected_aliases=aliases, normalized_rows=rows,
                                         metadata={"schemaVersion": 1, "genomeBuild": "GRCh37", "normalizationVersion": "synthetic",
                                                   "referenceSha256": digest, "sourceSha256": hashlib.sha256(self.source.read_bytes()).hexdigest(),
                                                   "excludedContigs": [], "observedAlleleCount": len(rows)})
            db.commit()
        return reader, reference, digest

    def test_ingestion_storage_staging_and_resolver_share_whole_call_semantics(self):
        reader, reference, digest = self._prepared_index()
        with sqlite3.connect(reader.agi_path) as db:
            records = list(db.execute("select rsid,genotype,observed_alleles,depth,genotype_quality from records"))
            metadata = {key: json.loads(value) for key, value in db.execute("select key,value from metadata")}
            self.assertEqual(metadata["array_call_summary"], {"total_rows": 15, "supported_single_base_calls": 11,
                                                           "incomplete_calls": 1, "unsupported_calls": 3})
            self.assertEqual(db.execute("select count(*) from normalized_alleles where pos in (30,40,50,60)").fetchone()[0], 0)
            self.assertEqual(len(records), 15)
            self.assertTrue(all(depth is None and quality is None for _, _, _, depth, quality in records))
        targets = [{"chrom": "1", "pos": pos, "ref": "A", "alt": "G" if pos == 10 else "C"}
                   for pos in (10, 20, 30, 40, 50, 60, 70)]
        resolution = reader.resolve_allele_targets(targets, reference_fasta=reference,
                                                   expected_reference_sha256=digest, deadline=time.monotonic() + 10)
        rows = resolution["results"]
        self.assertEqual([row["status"] for row in rows],
                         ["observed", "observed", "no_call", "unavailable", "unavailable", "unavailable", "reference_confirmed"])
        self.assertEqual([row["target_dosage"] for row in rows], [1, 2, None, None, None, None, 0])
        self.assertEqual(rows[0]["call_assessment"]["raw_genotype"], "A/G")
        self.assertEqual(rows[0]["call_assessment"]["alleles"], ["A", "G"])
        self.assertEqual(rows[0]["call_assessment"]["quality"], "unknown")
        self.assertEqual(rows[2]["source_records"][0]["genotype"], "./A")
        self.assertFalse(rows[3]["exact_match"])
        self.assertFalse(rows[4]["exact_match"])
        self.assertEqual(rows[6]["reference_scope"], "assayed_array_marker")

    def test_non_snp_public_targets_cannot_reuse_single_base_components(self):
        reader, reference, digest = self._prepared_index()
        targets = [{"chrom": "1", "pos": 10, "ref": ref, "alt": alt}
                   for ref, alt in (("AA", "G"), ("A", "GG"), ("A", "N"), ("A", "<DEL>"))]
        result = reader.resolve_allele_targets(targets, reference_fasta=reference,
                                              expected_reference_sha256=digest, deadline=time.monotonic() + 10)
        for row in result["results"]:
            self.assertEqual(row["status"], "unavailable")
            self.assertIsNone(row["target_dosage"])
            self.assertFalse(row["exact_match"])

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


if __name__ == "__main__":
    unittest.main()
