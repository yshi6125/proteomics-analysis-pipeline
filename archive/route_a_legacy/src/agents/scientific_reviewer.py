"""Bounded scientific quality-control review for generated pathway reports."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any, Callable

from .dataset_report_agent import LLMClient


MAX_REVIEW_ROUNDS = 2
VALID_ERROR_TYPES = {
    "dataset_fidelity", "unsupported_interpretation", "literature_fidelity",
    "citation_fidelity", "classification_mismatch", "completeness",
}


@dataclass(frozen=True)
class ReviewOutcome:
    final_report: str
    review_status: str
    review_rounds_used: int
    final_review: dict[str, Any]


class ScientificReviewerAgent:
    """Audit complete reports against supplied evidence without outside research."""

    def __init__(
        self,
        llm: LLMClient,
        prompt: str,
        *,
        provider: str,
        model: str,
        cache_dir: Path | None = None,
        cache_key: str | None = None,
        refresh_cache: bool = False,
    ):
        self.llm, self.prompt = llm, prompt
        self.provider, self.model = provider, model
        self.cache_dir, self.cache_key = cache_dir, cache_key
        self.refresh_cache = refresh_cache
        self.request_count = 0

    @staticmethod
    def _hash(value: Any) -> str:
        encoded = json.dumps(
            value, sort_keys=True, separators=(",", ":"), ensure_ascii=False
        ).encode("utf-8")
        return hashlib.sha256(encoded).hexdigest()

    def _cache_path(self, inputs: dict[str, Any]) -> Path | None:
        if self.cache_dir is None or self.cache_key is None:
            return None
        digest = self._hash(inputs)[:16]
        return self.cache_dir / f"scientific_review_{self.cache_key}_{digest}.json"

    def _validate(self, result: Any) -> dict[str, Any]:
        if not isinstance(result, dict) or result.get("status") not in {"pass", "revise"}:
            raise ValueError("Reviewer response status must be 'pass' or 'revise'.")
        errors = result.get("errors")
        if not isinstance(errors, list):
            raise ValueError("Reviewer response must contain an errors list.")
        if result["status"] == "pass" and errors:
            raise ValueError("Reviewer PASS response must not contain errors.")
        if result["status"] == "revise" and not errors:
            raise ValueError("Reviewer REVISE response must contain errors.")
        required = {"error_type", "cluster_id", "location", "claim_or_value", "reason", "required_fix"}
        for item in errors:
            if not isinstance(item, dict) or not required.issubset(item):
                raise ValueError("Each reviewer error must contain all required fields.")
            if item["error_type"] not in VALID_ERROR_TYPES:
                raise ValueError(f"Invalid reviewer error type: {item['error_type']}")
        return result

    def review(self, report: str, context: dict[str, Any]) -> dict[str, Any]:
        inputs = {
            "reviewer_prompt": self.prompt,
            "reviewer_provider": self.provider,
            "reviewer_model": self.model,
            "generated_report": report,
            **context,
        }
        path = self._cache_path(inputs)
        if path is not None and path.exists() and not self.refresh_cache:
            try:
                cached = json.loads(path.read_text(encoding="utf-8"))
                if cached.get("input_hash") == self._hash(inputs):
                    return self._validate(cached.get("result"))
            except (OSError, json.JSONDecodeError, AttributeError, ValueError):
                pass
        self.request_count += 1
        result = self._validate(self.llm.json_response(self.prompt, {
            "task": "review_report", "review_inputs": inputs,
        }))
        if path is not None:
            path.parent.mkdir(parents=True, exist_ok=True)
            temporary = path.with_suffix(path.suffix + ".tmp")
            temporary.write_text(json.dumps({
                "input_hash": self._hash(inputs), "result": result,
            }, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
            temporary.replace(path)
        return result

    def run(
        self,
        draft_report: str,
        context: dict[str, Any],
        revise: Callable[[str, list[dict[str, Any]]], str],
    ) -> ReviewOutcome:
        report = draft_report
        review: dict[str, Any] = {"status": "revise", "errors": []}
        for round_number in range(1, MAX_REVIEW_ROUNDS + 1):
            review = self.review(report, context)
            if review["status"] == "pass":
                return ReviewOutcome(report, "pass", round_number, review)
            if round_number < MAX_REVIEW_ROUNDS:
                report = revise(report, review["errors"])
        return ReviewOutcome(
            report, "manual_review_required", MAX_REVIEW_ROUNDS, review
        )
