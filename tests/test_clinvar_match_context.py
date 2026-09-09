from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from genomi.operations.registry import handlers_clinvar
from genomi.operations.registry.coerce import _with_context
from genomi.runtime import context as runtime_context
from genomi.runtime.paths import CLINVAR_MATCHES_NAME
from tests.support.runtime.genomi import GenomiRuntimeTestCase


class ClinvarMatchContextTests(GenomiRuntimeTestCase):
    def setUp(self):
        super().setUp()
        self.agi = self.genomi_home / "target" / "index.sqlite"
        self.agi.parent.mkdir(parents=True)
        self.agi.touch()
        self.run = {"agi_id": "target", "agi_path": str(self.agi), "matches": None}
        runtime_context.save_agi_to_registry(self.run)
        runtime_context.save_agi_to_registry({"agi_id": "other", "matches": "other.jsonl"})
        context = runtime_context.load_context()
        context["active_agi_id"] = "other"
        runtime_context.save_context(context)

    def test_completed_and_cached_matching_register_actual_output_without_switching_selection(self):
        for status in ("completed", "cached"):
            with self.subTest(status=status):
                output = self.agi.parent / f"custom-{status}.jsonl"
                output.touch()
                with (
                    patch.object(
                        handlers_clinvar,
                        "open_agi",
                        return_value=SimpleNamespace(agi_path=self.agi),
                    ),
                    patch.object(handlers_clinvar, "resolve_agi_record", return_value=self.run),
                    patch.object(
                        handlers_clinvar, "_with_context", return_value={"db": "evidence.sqlite"}
                    ),
                    patch.object(
                        handlers_clinvar.static_annotation,
                        "match_static_clinvar_from_active_genome_index",
                        return_value={"status": status, "output": str(output)},
                    ),
                ):
                    handlers_clinvar._clinvar_match({"agi_id": "target"})
                self.assertEqual(runtime_context.find_agi("target")["matches"], str(output))
                context = runtime_context.load_context()
                self.assertEqual(context["agis"]["target"]["matches"], str(output))
                self.assertEqual(context["active_agi_id"], "other")
                self.assertEqual(runtime_context.find_agi("other")["matches"], "other.jsonl")

    def test_failed_matching_does_not_register_an_artifact(self):
        with (
            patch.object(
                handlers_clinvar, "open_agi", return_value=SimpleNamespace(agi_path=self.agi)
            ),
            patch.object(handlers_clinvar, "resolve_agi_record", return_value=self.run),
            patch.object(handlers_clinvar, "_with_context", return_value={"db": "evidence.sqlite"}),
            patch.object(
                handlers_clinvar.static_annotation,
                "match_static_clinvar_from_active_genome_index",
                return_value={"status": "requires_library_install"},
            ),
        ):
            handlers_clinvar._clinvar_match({"agi_id": "target"})
        self.assertIsNone(runtime_context.find_agi("target")["matches"])

    def test_legacy_resolution_requires_approved_run_and_existing_artifact(self):
        artifact = self.agi.with_name(CLINVAR_MATCHES_NAME)
        with patch(
            "genomi.operations.registry.agi_access.resolve_agi_record", return_value=self.run
        ) as resolve:
            self.assertIsNone(_with_context({}, matches=True)["matches"])
            artifact.touch()
            self.assertEqual(_with_context({}, matches=True)["matches"], str(artifact))
            self.assertTrue(resolve.call_args.kwargs["require_approved"])
        with patch("genomi.operations.registry.agi_access.resolve_agi_record", return_value=None):
            self.assertNotIn("matches", _with_context({}, matches=True))
        self.assertIsNone(runtime_context.find_agi("target")["matches"])

    def test_recorded_and_explicit_paths_take_precedence_over_legacy_artifact(self):
        self.agi.with_name(CLINVAR_MATCHES_NAME).touch()
        with patch(
            "genomi.operations.registry.agi_access.resolve_agi_record",
            return_value={**self.run, "matches": "recorded.jsonl"},
        ):
            self.assertEqual(_with_context({}, matches=True)["matches"], "recorded.jsonl")
            self.assertEqual(
                _with_context({"matches": "explicit.jsonl"}, matches=True)["matches"],
                "explicit.jsonl",
            )
