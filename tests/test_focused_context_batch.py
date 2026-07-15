from __future__ import annotations

import unittest
from unittest.mock import patch

from genomi.operations import OperationError
from genomi.operations.registry.handlers_evidence_phenotype import (
    _evidence_gather_focused_context,
)


class FocusedContextBatchTests(unittest.TestCase):
    def test_gathers_mixed_targets_in_one_bounded_operation(self) -> None:
        with patch("genomi.operations.registry.table.call_operation") as call:
            call.side_effect = [
                {"resolved": "rs123"},
                {"query": {"gene": "GENE1"}},
                OperationError("evidence_unavailable", "Medication evidence is unavailable."),
            ]
            result = _evidence_gather_focused_context(
                {
                    "rsids": ["rs123"],
                    "genes": ["GENE1"],
                    "medications": ["example-drug"],
                }
            )

        self.assertEqual(result["target_count"], 3)
        self.assertEqual(
            [item["status"] for item in result["results"]],
            ["completed", "completed", "unavailable"],
        )
        self.assertEqual(call.call_count, 3)
        self.assertEqual(call.call_args_list[0].args[0], "variant.resolve")
        self.assertEqual(call.call_args_list[1].args[0], "variant.gather_gene_context")
        self.assertEqual(call.call_args_list[2].args[0], "pharmacogenomics.review_medication")

    def test_requires_a_target_and_caps_combined_targets(self) -> None:
        with self.assertRaises(OperationError) as missing:
            _evidence_gather_focused_context({})
        self.assertEqual(missing.exception.code, "needs_input")

        with self.assertRaises(OperationError) as too_many:
            _evidence_gather_focused_context({"rsids": [f"rs{index}" for index in range(9)]})
        self.assertEqual(too_many.exception.code, "invalid_params")


if __name__ == "__main__":
    unittest.main()
