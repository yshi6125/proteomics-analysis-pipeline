from datetime import date
import json
from urllib import error, parse

import pytest

from src.agents.dataset_report_agent import (
    DatasetReportAgent, GeminiClient, OpenAIResponsesClient, PubMedClient, Study,
    dataset_evidence_strength, disease_relevance, load_pathway_evidence,
)
from src import generate_pathway_report


def cluster(l1=0, l2=0, l3=0, l4=0, l5=0):
    return {"level_counts": {"1": l1, "2": l2, "3": l3, "4": l4, "5": l5}}


@pytest.mark.parametrize(("value", "expected"), [
    (cluster(l1=2), "STRONG"), (cluster(l1=1, l2=1, l3=1), "STRONG"),
    (cluster(l1=1, l2=1), "MODERATE"), (cluster(l2=2), "MODERATE"),
    (cluster(l3=1, l4=10), "EXPLORATORY"),
])
def test_dataset_strength(value, expected):
    assert dataset_evidence_strength(value) == expected


def test_disease_relevance_is_literature_only():
    direct = [
        {"is_primary_study": True, "direct_disease_support": True,
         "evidence_type": kind} for kind in ("mechanistic", "clinical", "associative")
    ]
    assert disease_relevance(direct) == "HIGH"
    assert disease_relevance(direct[:2]) == "MODERATE"
    assert disease_relevance([]) == "LOW"


def test_pubmed_query_is_exact_five_year_window():
    query = PubMedClient.build_query("fatty acid oxidation", "type 2 diabetes", date(2026, 8, 14))
    assert '"2021/08/14"[Date - Publication] : "2026/08/14"[Date - Publication]' in query
    assert '"fatty acid oxidation"[Title/Abstract]' in query


class FakeClock:
    def __init__(self):
        self.now = 0.0
        self.sleeps = []

    def monotonic(self):
        return self.now

    def sleep(self, seconds):
        self.sleeps.append(seconds)
        self.now += seconds


class HTTPResponse:
    def __init__(self, body=b"ok"):
        self.body = body

    def __enter__(self):
        return self

    def __exit__(self, *_args):
        return None

    def read(self):
        return self.body


def http_error(code, retry_after=None):
    headers = {} if retry_after is None else {"Retry-After": str(retry_after)}
    return error.HTTPError("https://ncbi.test", code, "failure", headers, None)


def test_pubmed_normal_requests_are_throttled_without_api_key(monkeypatch):
    monkeypatch.delenv("NCBI_API_KEY", raising=False)
    clock = FakeClock()
    client = PubMedClient(
        opener=lambda *_args, **_kwargs: HTTPResponse(),
        sleep=clock.sleep, monotonic=clock.monotonic,
    )
    client._get("esearch.fcgi", {"db": "pubmed"})
    client._get("efetch.fcgi", {"db": "pubmed"})
    assert clock.sleeps == pytest.approx([0.4])


def test_pubmed_429_retries_then_succeeds(caplog, monkeypatch):
    monkeypatch.delenv("NCBI_API_KEY", raising=False)
    outcomes = [http_error(429), HTTPResponse()]
    clock = FakeClock()

    def opener(*_args, **_kwargs):
        outcome = outcomes.pop(0)
        if isinstance(outcome, Exception):
            raise outcome
        return outcome

    client = PubMedClient(
        opener=opener, sleep=clock.sleep, monotonic=clock.monotonic,
    )
    assert client._get("esearch.fcgi", {}) == b"ok"
    assert clock.sleeps == [2.0]
    assert "PubMed rate limited (429); retry 1/4 in 2 seconds." in caplog.text


def test_pubmed_retry_after_overrides_shorter_backoff(monkeypatch):
    monkeypatch.delenv("NCBI_API_KEY", raising=False)
    outcomes = [http_error(429, retry_after=7), HTTPResponse()]
    clock = FakeClock()

    def opener(*_args, **_kwargs):
        outcome = outcomes.pop(0)
        if isinstance(outcome, Exception):
            raise outcome
        return outcome

    client = PubMedClient(
        opener=opener, sleep=clock.sleep, monotonic=clock.monotonic,
    )
    client._get("esearch.fcgi", {})
    assert clock.sleeps == [7.0]


def test_pubmed_transient_retries_exhausted(monkeypatch):
    monkeypatch.delenv("NCBI_API_KEY", raising=False)
    calls = 0
    clock = FakeClock()

    def opener(*_args, **_kwargs):
        nonlocal calls
        calls += 1
        raise http_error(503)

    client = PubMedClient(
        opener=opener, sleep=clock.sleep, monotonic=clock.monotonic,
    )
    with pytest.raises(RuntimeError, match="after 4 retries"):
        client._get("efetch.fcgi", {})
    assert calls == 5
    assert clock.sleeps == [2.0, 4.0, 8.0, 16.0]


@pytest.mark.parametrize("status_code", [400, 401, 403, 404])
def test_pubmed_permanent_http_error_fails_immediately(status_code, monkeypatch):
    monkeypatch.delenv("NCBI_API_KEY", raising=False)
    calls = 0
    clock = FakeClock()

    def opener(*_args, **_kwargs):
        nonlocal calls
        calls += 1
        raise http_error(status_code)

    client = PubMedClient(
        opener=opener, sleep=clock.sleep, monotonic=clock.monotonic,
    )
    with pytest.raises(RuntimeError, match=f"HTTP {status_code}"):
        client._get("esearch.fcgi", {})
    assert calls == 1
    assert clock.sleeps == []


def test_pubmed_api_key_enables_higher_rate_and_is_sent(monkeypatch):
    monkeypatch.setenv("NCBI_API_KEY", "secret-key")
    clock = FakeClock()
    urls = []

    def opener(url, **_kwargs):
        urls.append(url)
        return HTTPResponse()

    client = PubMedClient(
        opener=opener, sleep=clock.sleep, monotonic=clock.monotonic,
    )
    client._get("esearch.fcgi", {})
    client._get("efetch.fcgi", {})
    assert clock.sleeps == pytest.approx([0.1])
    assert all(
        parse.parse_qs(parse.urlparse(url).query)["api_key"] == ["secret-key"]
        for url in urls
    )


@pytest.mark.parametrize("publication_type", [
    "Review", "Systematic Review", "Editorial", "Comment", "Perspective",
])
def test_background_publication_types_are_not_primary(publication_type):
    study = Study(
        pmid="1", title="title", abstract="", journal="journal", year=2025,
        authors=(), publication_types=(publication_type,),
    )
    assert study.is_non_primary_publication


def test_loader_rejects_missing_clusters(tmp_path):
    path = tmp_path / "evidence.json"
    path.write_text(json.dumps({"number_of_clusters": 1, "clusters": []}))
    with pytest.raises(ValueError, match="does not match"):
        load_pathway_evidence(path)


def test_preview_builds_queries_without_pubmed_retrieval_or_report_generation():
    class ThemeOnlyLLM:
        calls = []

        def json_response(self, _instructions, payload):
            self.calls.append(payload["task"])
            assert payload["task"] == "name_biological_themes"
            return {"themes": [{"cluster_id": "C001", "biological_theme": "mitochondrial respiration"}]}

    class NoRetrievalPubMed(PubMedClient):
        def search(self, *args, **kwargs):
            raise AssertionError("Preview must not call PubMed")

    evidence = {"clusters": [{
        "cluster_id": "C001",
        "representative_pathway": {"name": "Cellular Respiration", "library": "GO"},
        "member_pathways": [],
        "genes": [{"gene": "NDUFA10", "level": 2}],
    }]}
    llm = ThemeOnlyLLM()
    agent = DatasetReportAgent(llm, NoRetrievalPubMed(), "prompt")

    preview = agent.preview_queries(evidence, "Type 2 diabetes")

    assert llm.calls == ["name_biological_themes"]
    assert preview[0]["cluster_id"] == "C001"
    assert preview[0]["representative_pathway"] == "Cellular Respiration"
    assert preview[0]["biological_search_theme"] == "mitochondrial respiration"
    assert '"Type 2 diabetes"[Title/Abstract]' in preview[0]["pubmed_query"]


def test_gemini_requires_api_key(monkeypatch):
    monkeypatch.delenv("GEMINI_API_KEY", raising=False)
    with pytest.raises(ValueError, match="GEMINI_API_KEY"):
        GeminiClient("gemini-test")


def test_gemini_response_parsing(monkeypatch):
    class Models:
        def generate_content(self, *, model, contents):
            assert model == "gemini-test"
            assert "Input JSON:" in contents
            return type("Response", (), {"text": '```json\n{"themes": []}\n```'})()

    client = type("Client", (), {"models": Models()})()
    gemini = GeminiClient("gemini-test", api_key="test-key", client=client)
    assert gemini.json_response("instructions", {"task": "test"}) == {"themes": []}


class GeminiAPIError(Exception):
    def __init__(self, code):
        self.code = code
        super().__init__(f"Gemini API error {code}")


def test_gemini_retries_transient_failure_then_succeeds(caplog):
    class Models:
        calls = 0

        def generate_content(self, **_kwargs):
            self.calls += 1
            if self.calls == 1:
                raise GeminiAPIError(503)
            return type("Response", (), {"text": "success"})()

    models = Models()
    sleeps = []
    client = type("Client", (), {"models": models})()
    gemini = GeminiClient(
        "gemini-test", api_key="test-key", client=client, sleep=sleeps.append
    )

    assert gemini.generate("prompt") == "success"
    assert models.calls == 2
    assert sleeps == [2]
    assert "Gemini request failed with 503; retry 1/4 in 2 seconds." in caplog.text


def test_gemini_maximum_retries_exhausted():
    original = GeminiAPIError(504)

    class Models:
        calls = 0

        def generate_content(self, **_kwargs):
            self.calls += 1
            raise original

    models = Models()
    sleeps = []
    client = type("Client", (), {"models": models})()
    gemini = GeminiClient(
        "gemini-test", api_key="test-key", client=client, sleep=sleeps.append
    )

    with pytest.raises(GeminiAPIError) as raised:
        gemini.generate("prompt")
    assert raised.value is original
    assert models.calls == 5
    assert sleeps == [2, 4, 8, 16]
    assert any("failed after 4 retries" in note for note in original.__notes__)


@pytest.mark.parametrize("status_code", [400, 401, 403])
def test_gemini_non_retryable_error_fails_immediately(status_code):
    class Models:
        calls = 0

        def generate_content(self, **_kwargs):
            self.calls += 1
            raise GeminiAPIError(status_code)

    models = Models()
    sleeps = []
    client = type("Client", (), {"models": models})()
    gemini = GeminiClient(
        "gemini-test", api_key="test-key", client=client, sleep=sleeps.append
    )

    with pytest.raises(GeminiAPIError):
        gemini.generate("prompt")
    assert models.calls == 1
    assert sleeps == []


def test_gemini_daily_quota_exhaustion_is_not_retried():
    error = GeminiAPIError(429)
    error.args = ("429 daily quota exhausted",)

    class Models:
        calls = 0

        def generate_content(self, **_kwargs):
            self.calls += 1
            raise error

    models = Models()
    sleeps = []
    gemini = GeminiClient(
        "gemini-test", api_key="test-key",
        client=type("Client", (), {"models": models})(), sleep=sleeps.append,
    )
    with pytest.raises(GeminiAPIError):
        gemini.generate("prompt")
    assert models.calls == 1
    assert sleeps == []
    assert any("quota exhaustion" in note for note in error.__notes__)


@pytest.mark.parametrize("provider", ["openai", "gemini"])
def test_provider_selection(monkeypatch, provider):
    sentinel = object()
    client_name = "OpenAIResponsesClient" if provider == "openai" else "GeminiClient"
    monkeypatch.setattr(generate_pathway_report, client_name, lambda model: (sentinel, model))
    assert generate_pathway_report.build_llm_client(provider, "model-name") == (
        sentinel, "model-name"
    )


def test_existing_openai_output_text_parsing_is_unchanged(monkeypatch):
    class Response:
        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return None

        def read(self):
            return b'{"output_text":"{\\"ok\\": true}"}'

    monkeypatch.setattr(
        "src.agents.dataset_report_agent.request.urlopen",
        lambda *_args, **_kwargs: Response(),
    )
    client = OpenAIResponsesClient("openai-test", api_key="test-key")
    assert client.json_response("instructions", {"task": "test"}) == {"ok": True}


def report_evidence(number=5):
    return {
        "metadata": {"comparison": "A_vs_B"},
        "clusters": [
            {
                "cluster_id": f"C{index:03d}",
                "direction": "decreased",
                "representative_pathway": {"name": f"Pathway {index}", "library": "GO"},
                "member_pathways": [{"name": f"Pathway {index}", "library": "GO"}],
                "genes": [{"gene": f"GENE{index}", "level": 1}],
                "level_counts": {"1": 1, "2": 0, "3": 1, "4": 0, "5": 0},
                "best_adjusted_p": 0.01,
                "representative_gene_ratio": 0.2,
            }
            for index in range(1, number + 1)
        ],
    }


class BatchPubMed(PubMedClient):
    def search(self, theme, disease, **_kwargs):
        cluster_number = theme.removeprefix("theme ")
        study = Study(
            pmid=f"PMID{cluster_number}", title=f"Study {cluster_number}",
            abstract="abstract", journal="journal", year=2025,
            authors=("Author",), publication_types=("Journal Article",),
        )
        return self.build_query(theme, disease, date(2026, 8, 14)), [study]


class BatchLLM:
    def __init__(self):
        self.payloads = []

    def json_response(self, _prompt, payload):
        self.payloads.append(payload)
        if payload["task"] == "name_biological_themes":
            return {"themes": [
                {"cluster_id": item["cluster_id"],
                 "biological_theme": f"theme {int(item['cluster_id'][1:])}"}
                for item in payload["clusters"]
            ]}
        if payload["task"] == "assess_literature":
            return {"cluster_assessments": [
                {
                    "cluster_id": item["cluster_id"],
                    "study_annotations": [{
                        "pmid": record["pmid"], "is_primary_study": True,
                        "direct_disease_support": True,
                        "closely_related_model_support": False,
                        "evidence_type": "clinical", "tissue_or_model": "human",
                        "finding": "association",
                    } for record in item["pubmed_records"]],
                    "major_findings": [], "agreement": "consistent",
                    "mechanisms": [], "common_tissues_or_models": ["human"],
                }
                for item in payload["clusters"]
            ]}
        assert payload["task"] == "write_report"
        return {
            "executive_summary": "summary",
            "cluster_sections": [
                {"cluster_id": item["cluster"]["cluster_id"], "markdown": "section"}
                for item in payload["clusters"]
            ],
            "final_biological_interpretation": "interpretation",
        }


def test_literature_is_batched_without_cluster_or_record_leakage():
    llm = BatchLLM()
    agent = DatasetReportAgent(
        llm, BatchPubMed(), "prompt", literature_batch_size=2
    )
    evidence = report_evidence(5)

    agent.generate(evidence, "Disease")

    tasks = [payload["task"] for payload in llm.payloads]
    assert tasks == [
        "name_biological_themes", "assess_literature", "assess_literature",
        "assess_literature", "write_report",
    ]
    batches = [p["clusters"] for p in llm.payloads if p["task"] == "assess_literature"]
    assert [[item["cluster_id"] for item in batch] for batch in batches] == [
        ["C001", "C002"], ["C003", "C004"], ["C005"],
    ]
    for batch in batches:
        for item in batch:
            expected_pmid = f"PMID{int(item['cluster_id'][1:])}"
            assert [record["pmid"] for record in item["pubmed_records"]] == [expected_pmid]
    assert agent.request_counts == {
        "biological_themes": 1, "literature_assessment": 3, "report_writing": 1,
        "revision": 0,
    }


def test_cache_avoids_intermediate_calls_and_refresh_bypasses_it(tmp_path):
    evidence = report_evidence(2)
    first_llm = BatchLLM()
    first = DatasetReportAgent(
        first_llm, BatchPubMed(), "prompt", literature_batch_size=2,
        cache_dir=tmp_path, cache_key="comparison_disease",
    )
    first.generate(evidence, "Disease")
    assert first.request_counts == {
        "biological_themes": 1, "literature_assessment": 1, "report_writing": 1,
        "revision": 0,
    }

    cached_llm = BatchLLM()
    cached = DatasetReportAgent(
        cached_llm, BatchPubMed(), "prompt", literature_batch_size=2,
        cache_dir=tmp_path, cache_key="comparison_disease",
    )
    cached.generate(evidence, "Disease")
    assert [payload["task"] for payload in cached_llm.payloads] == ["write_report"]
    assert cached.request_counts == {
        "biological_themes": 0, "literature_assessment": 0, "report_writing": 1,
        "revision": 0,
    }

    refreshed_llm = BatchLLM()
    refreshed = DatasetReportAgent(
        refreshed_llm, BatchPubMed(), "prompt", literature_batch_size=2,
        cache_dir=tmp_path, cache_key="comparison_disease", refresh_cache=True,
    )
    refreshed.generate(evidence, "Disease")
    assert refreshed.request_counts == {
        "biological_themes": 1, "literature_assessment": 1, "report_writing": 1,
        "revision": 0,
    }


def test_changed_pathway_evidence_invalidates_theme_cache(tmp_path):
    evidence = report_evidence(2)
    first = DatasetReportAgent(
        BatchLLM(), BatchPubMed(), "prompt", literature_batch_size=2,
        cache_dir=tmp_path, cache_key="comparison_disease",
    )
    first.generate(evidence, "Disease")

    changed = json.loads(json.dumps(evidence))
    changed["clusters"][0]["genes"].append({"gene": "NEW_GENE", "level": 4})
    changed_llm = BatchLLM()
    changed_agent = DatasetReportAgent(
        changed_llm, BatchPubMed(), "prompt", literature_batch_size=2,
        cache_dir=tmp_path, cache_key="comparison_disease",
    )
    changed_agent.generate(changed, "Disease")

    assert changed_agent.request_counts["biological_themes"] == 1
    # Literature assessment may remain reusable when the newly generated theme,
    # PubMed query, and retrieved records are byte-for-byte unchanged.
    assert changed_agent.request_counts["literature_assessment"] == 0
