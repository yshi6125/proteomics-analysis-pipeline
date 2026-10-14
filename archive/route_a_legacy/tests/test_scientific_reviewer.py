import pytest

from src.agents.scientific_reviewer import (
    MAX_REVIEW_ROUNDS, ScientificReviewerAgent,
)


def reviewer_error(error_type="dataset_fidelity", claim="wrong value"):
    return {
        "error_type": error_type,
        "cluster_id": "C001",
        "location": "C001 dataset evidence",
        "claim_or_value": claim,
        "reason": "conflicts with supplied evidence",
        "required_fix": "restore the supplied value",
    }


class SequencedReviewerLLM:
    def __init__(self, results):
        self.results = list(results)
        self.payloads = []

    def json_response(self, _prompt, payload):
        self.payloads.append(payload)
        return self.results.pop(0)


def context():
    return {
        "pathway_evidence": {"clusters": [
            {"cluster_id": "C001"}, {"cluster_id": "C002"},
        ]},
        "biological_themes": {"C001": "respiration", "C002": "translation"},
        "retrieved_pubmed": [
            {"cluster_id": "C001", "pubmed_records": [{"pmid": "1"}]},
            {"cluster_id": "C002", "pubmed_records": [{"pmid": "2"}]},
        ],
        "literature_assessments": {
            "C001": {"study_annotations": [{"pmid": "1"}]},
            "C002": {"study_annotations": [{"pmid": "2"}]},
        },
        "deterministic_clusters": [
            {"cluster": {"cluster_id": "C001"}, "dataset_evidence_strength": "MODERATE"},
            {"cluster": {"cluster_id": "C002"}, "dataset_evidence_strength": "STRONG"},
        ],
    }


def make_reviewer(llm, tmp_path=None, refresh=False):
    return ScientificReviewerAgent(
        llm, "reviewer prompt", provider="gemini", model="reviewer-model",
        cache_dir=tmp_path, cache_key="A_vs_B_disease", refresh_cache=refresh,
    )


def test_clean_report_passes_without_revision_and_receives_all_clusters():
    llm = SequencedReviewerLLM([{"status": "pass", "errors": []}])
    reviewer = make_reviewer(llm)
    revisions = []
    outcome = reviewer.run(
        "clean report", context(),
        lambda report, errors: revisions.append((report, errors)) or report,
    )
    assert outcome.review_status == "pass"
    assert outcome.review_rounds_used == 1
    assert revisions == []
    assert reviewer.request_count == 1
    supplied = llm.payloads[0]["review_inputs"]
    assert [item["cluster_id"] for item in supplied["pathway_evidence"]["clusters"]] == [
        "C001", "C002",
    ]


@pytest.mark.parametrize(("error_type", "claim"), [
    ("dataset_fidelity", "Total proteins: 99"),
    ("classification_mismatch", "Dataset Evidence Strength: STRONG"),
    ("unsupported_interpretation", "pathway activated"),
    ("citation_fidelity", "PMID: 999999"),
    ("citation_fidelity", "PMID 2 listed under C001"),
    ("completeness", "C002 missing"),
])
def test_material_report_errors_request_revision(error_type, claim):
    llm = SequencedReviewerLLM([{
        "status": "revise", "errors": [reviewer_error(error_type, claim)],
    }])
    result = make_reviewer(llm).review(f"report with {claim}", context())
    assert result["status"] == "revise"
    assert result["errors"][0]["error_type"] == error_type


def test_successful_revision_is_reviewed_again_and_passes():
    issue = reviewer_error()
    llm = SequencedReviewerLLM([
        {"status": "revise", "errors": [issue]},
        {"status": "pass", "errors": []},
    ])
    revision_calls = []

    def revise(report, errors):
        revision_calls.append((report, errors))
        return "corrected report"

    reviewer = make_reviewer(llm)
    outcome = reviewer.run("draft", context(), revise)
    assert outcome.final_report == "corrected report"
    assert outcome.review_status == "pass"
    assert outcome.review_rounds_used == 2
    assert reviewer.request_count == 2
    assert len(revision_calls) == 1


def test_unresolved_second_review_requires_manual_review_and_stops():
    issue = reviewer_error()
    llm = SequencedReviewerLLM([
        {"status": "revise", "errors": [issue]},
        {"status": "revise", "errors": [issue]},
        {"status": "pass", "errors": []},
    ])
    revision_calls = []
    reviewer = make_reviewer(llm)
    outcome = reviewer.run(
        "draft", context(),
        lambda *_args: revision_calls.append(True) or "still wrong",
    )
    assert MAX_REVIEW_ROUNDS == 2
    assert outcome.review_status == "manual_review_required"
    assert outcome.final_report == "still wrong"
    assert reviewer.request_count == 2
    assert len(revision_calls) == 1
    assert len(llm.results) == 1


def test_reviewer_cache_reuse_report_invalidation_and_refresh(tmp_path):
    passed = {"status": "pass", "errors": []}
    first_llm = SequencedReviewerLLM([passed])
    assert make_reviewer(first_llm, tmp_path).review("report one", context()) == passed

    cached_llm = SequencedReviewerLLM([])
    cached = make_reviewer(cached_llm, tmp_path)
    assert cached.review("report one", context()) == passed
    assert cached.request_count == 0

    changed_llm = SequencedReviewerLLM([passed])
    changed = make_reviewer(changed_llm, tmp_path)
    changed.review("report two", context())
    assert changed.request_count == 1

    refresh_llm = SequencedReviewerLLM([passed])
    refreshed = make_reviewer(refresh_llm, tmp_path, refresh=True)
    refreshed.review("report one", context())
    assert refreshed.request_count == 1
