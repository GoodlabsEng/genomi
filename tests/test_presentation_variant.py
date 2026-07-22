from __future__ import annotations

import unittest

from genomi.interfaces.presentation import present_result


class VariantPresentationTests(unittest.TestCase):
    def test_variant_resolution_preserves_interpretation_evidence(self) -> None:
        raw_variant = {
            "query": {
                "text": "rs762551",
                "genome_build": "GRCh38",
                "agi_id": "agi-private-id",
                "include_active_genome_index": True,
                "limit": 20,
            },
            "resolved_targets": [{"target_type": "rsid", "rsid": "rs762551"}],
            "sample_context": {
                "searched_active_genome_indexes": [
                    {
                        "agi_id": "agi-private-id",
                        "sample_slug": "sample-private-id",
                        "availability": {"agi_path": "/tmp/private/sample.sqlite"},
                    }
                ],
                "searched_known_active_genome_indexes": False,
                "count": 1,
                "matches": [
                    {
                        "target": "rs762551",
                        "rsid": "rs762551",
                        "chrom": "15",
                        "pos": 75041917,
                        "ref": "A",
                        "alt": "C",
                        "genotype": "AC",
                        "observed_alleles": ["A", "C"],
                        "filter": "PASS",
                        "is_variant": True,
                        "record_kind": "variant",
                        "depth": 28,
                        "genotype_quality": 99,
                        "agi_source_format": "23andme",
                        "agi_source_kind": "genotype_array",
                        "selection": "active_genome_index",
                        "agi_id": "agi-private-id",
                        "sample_slug": "sample-private-id",
                        "sample_name": "Example Person",
                        "sample_raw": "AC",
                        "offset": 42,
                    }
                ],
            },
            "public_context": {
                "evidence_stores": [
                    {
                        "label": "agi:agi-private-id",
                        "available": True,
                        "size_bytes": 1234,
                    }
                ],
                "clinvar_by_rsid": [
                    {
                        "rsid": "rs762551",
                        "clinical_significance": "drug response",
                        "raw_json": {"internal": "large"},
                        "source_path": "/tmp/private/clinvar.sqlite",
                    }
                ],
                "clinvar_by_allele": [],
                "clinvar_by_locus": [],
                "population_frequencies": [],
                "reviewed_research": [],
            },
            "support_context": {
                "genotype_support": [
                    {
                        "agi_id": "agi-private-id",
                        "sample_slug": "sample-private-id",
                        "chrom": "15",
                        "pos": 75041917,
                        "ref": "A",
                        "alt": "C",
                        "genome_build": "GRCh38",
                        "support_status": "supported",
                        "evidence_class": "genotype_support_supported",
                        "genotype": "AC",
                        "sample_observation": {
                            "target_alt_observed": True,
                            "alt_allele_count": 1,
                        },
                        "site_observation": {
                            "reference_call_supported": False,
                            "site_status": "variant_observed",
                        },
                    }
                ]
            },
            "target_inventory": {
                "target_count": 1,
                "has_sample_matches": True,
                "has_genotype_support": True,
            },
            "unanswered_answer_components": [],
            "warnings": [],
        }

        presented = present_result("variant.resolve", raw_variant)

        self.assertEqual(
            presented["query"],
            {"text": "rs762551", "genome_build": "GRCh38"},
        )
        self.assertEqual(presented["sample_context"]["active_genome_indexes_searched"], 1)
        self.assertEqual(
            presented["sample_context"]["matches"],
            [
                {
                    "target": "rs762551",
                    "rsid": "rs762551",
                    "chrom": "15",
                    "pos": 75041917,
                    "ref": "A",
                    "alt": "C",
                    "genotype": "AC",
                    "observed_alleles": ["A", "C"],
                    "filter": "PASS",
                    "is_variant": True,
                    "record_kind": "variant",
                    "depth": 28,
                    "genotype_quality": 99,
                    "agi_source_format": "23andme",
                    "agi_source_kind": "genotype_array",
                    "selection": "active_genome_index",
                }
            ],
        )
        self.assertEqual(
            presented["support_context"]["genotype_support"][0]["sample_observation"],
            {"target_alt_observed": True, "alt_allele_count": 1},
        )
        self.assertEqual(
            presented["public_context"]["clinvar_by_rsid"],
            [{"rsid": "rs762551", "clinical_significance": "drug response"}],
        )

    def test_focused_context_composes_normal_presenters(self) -> None:
        presented = present_result(
            "variant.gather_focused_context",
            {
                "status": "completed",
                "target_count": 2,
                "results": [
                    {
                        "target_type": "rsid",
                        "target": "rs762551",
                        "operation": "variant.resolve",
                        "status": "completed",
                        "result": {
                            "query": {"text": "rs762551", "genome_build": "GRCh38"},
                            "sample_context": {
                                "count": 1,
                                "searched_active_genome_indexes": [{}],
                                "matches": [{"rsid": "rs762551", "genotype": "AC"}],
                            },
                            "public_context": {},
                            "support_context": {},
                        },
                    },
                    {
                        "target_type": "medication",
                        "target": "example-drug",
                        "operation": "pharmacogenomics.review_medication",
                        "status": "unavailable",
                        "error": {
                            "code": "evidence_unavailable",
                            "message": "Medication evidence is unavailable.",
                        },
                    },
                ],
            },
        )

        self.assertEqual(presented["target_count"], 2)
        self.assertEqual(
            presented["results"][0]["result"]["sample_context"]["matches"],
            [{"rsid": "rs762551", "genotype": "AC"}],
        )
        self.assertEqual(
            presented["results"][1]["error"]["code"],
            "evidence_unavailable",
        )


if __name__ == "__main__":
    unittest.main()
