"""Generate a disease-focused report from frozen pathway-evidence JSON.

Dataset strength and disease relevance are assigned by explicit Python rules.  The
LLM is used only to name themes, assess retrieved papers, and write cautious prose.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import date, datetime, timezone
from email.utils import parsedate_to_datetime
import hashlib
import json
import logging
import os
from pathlib import Path
import re
import time
import threading
from typing import Any, Callable, Iterable, Protocol
from urllib import error, parse, request
import xml.etree.ElementTree as ET


PUBMED_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
LOGGER = logging.getLogger(__name__)
GEMINI_TRANSIENT_STATUS_CODES = frozenset({429, 500, 502, 503, 504})
GEMINI_RETRY_DELAYS = (2, 4, 8, 16)
DEFAULT_ABSTRACT_MAX_CHARS = 6000
NCBI_RETRY_DELAYS = (2, 4, 8, 16)
NCBI_TRANSIENT_STATUS_CODES = frozenset({429, 500, 502, 503, 504})
NCBI_DELAY_WITHOUT_KEY = 0.4
NCBI_DELAY_WITH_KEY = 0.1
REQUIRED_CLUSTER_FIELDS = {
    "cluster_id", "direction", "representative_pathway", "member_pathways",
    "genes", "level_counts", "best_adjusted_p", "representative_gene_ratio",
}


class LLMClient(Protocol):
    """Common structured-output interface consumed by DatasetReportAgent."""

    def json_response(
        self, instructions: str, payload: dict[str, Any]
    ) -> dict[str, Any]:
        """Return one JSON object for the supplied instructions and payload."""
        ...


@dataclass(frozen=True)
class Study:
    pmid: str
    title: str
    abstract: str
    journal: str
    year: int | None
    authors: tuple[str, ...]
    publication_types: tuple[str, ...]

    @property
    def is_non_primary_publication(self) -> bool:
        """Identify publication types that cannot count as primary research."""
        excluded = ("review", "editorial", "comment", "perspective")
        return any(
            any(label in item.lower() for label in excluded)
            for item in self.publication_types
        )

    @property
    def is_review(self) -> bool:
        """Retained for explicit review metadata in the model input."""
        return any("review" in item.lower() for item in self.publication_types)


def load_pathway_evidence(path: Path) -> dict[str, Any]:
    """Load and validate the stable evidence-builder contract."""
    try:
        evidence = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        raise FileNotFoundError(f"Pathway evidence not found: {path}") from None
    if not isinstance(evidence, dict) or not isinstance(evidence.get("clusters"), list):
        raise ValueError("Evidence must be an object containing a clusters list.")
    clusters = evidence["clusters"]
    if evidence.get("number_of_clusters") != len(clusters):
        raise ValueError("number_of_clusters does not match the clusters list.")
    ids: list[str] = []
    for index, cluster in enumerate(clusters):
        if not isinstance(cluster, dict):
            raise ValueError(f"Cluster {index} must be an object.")
        missing = sorted(REQUIRED_CLUSTER_FIELDS.difference(cluster))
        if missing:
            raise ValueError(f"Cluster {index} is missing fields: {missing}")
        ids.append(str(cluster["cluster_id"]))
        counts = cluster["level_counts"]
        if not isinstance(counts, dict) or any(str(level) not in counts for level in range(1, 6)):
            raise ValueError(f"Cluster {ids[-1]} must contain Level 1-5 counts.")
    if len(ids) != len(set(ids)):
        raise ValueError("Cluster IDs must be unique.")
    return evidence


def dataset_evidence_strength(cluster: dict[str, Any]) -> str:
    """Apply the frozen, deterministic dataset-strength rules."""
    counts = {level: int(cluster["level_counts"].get(str(level), 0)) for level in range(1, 6)}
    level_1 = counts[1]
    additional_2_3 = counts[2] + counts[3]
    levels_1_3 = level_1 + additional_2_3
    if level_1 >= 2 or (level_1 >= 1 and additional_2_3 >= 2):
        return "STRONG"
    if levels_1_3 >= 2:
        return "MODERATE"
    return "EXPLORATORY"


def disease_relevance(annotations: Iterable[dict[str, Any]]) -> str:
    """Classify literature evidence independently of dataset evidence."""
    primary = [a for a in annotations if a.get("is_primary_study")]
    direct = [a for a in primary if a.get("direct_disease_support")]
    substantive = [
        a for a in direct
        if a.get("evidence_type") in {"mechanistic", "clinical", "disease_progression"}
    ]
    if len(direct) >= 3 and len(substantive) >= 2:
        return "HIGH"
    related = any(a.get("closely_related_model_support") for a in primary)
    if direct or related:
        return "MODERATE"
    return "LOW"


class PubMedClient:
    """Small dependency-free PubMed E-utilities client."""

    def __init__(
        self,
        *,
        email: str | None = None,
        api_key: str | None = None,
        timeout: float = 30.0,
        opener: Callable[..., Any] = request.urlopen,
        sleep: Callable[[float], None] = time.sleep,
        monotonic: Callable[[], float] = time.monotonic,
    ):
        self.email = email or os.getenv("NCBI_EMAIL")
        self.api_key = api_key or os.getenv("NCBI_API_KEY")
        self.timeout = timeout
        self._opener = opener
        self._sleep = sleep
        self._monotonic = monotonic
        self._minimum_interval = (
            NCBI_DELAY_WITH_KEY if self.api_key else NCBI_DELAY_WITHOUT_KEY
        )
        self._last_request_at: float | None = None
        self._request_lock = threading.Lock()

    def _throttle(self) -> None:
        """Space every E-utilities request according to NCBI rate guidance."""
        with self._request_lock:
            now = self._monotonic()
            if self._last_request_at is not None:
                remaining = self._minimum_interval - (now - self._last_request_at)
                if remaining > 0:
                    self._sleep(remaining)
                    now = self._monotonic()
            self._last_request_at = now

    @staticmethod
    def _retry_after_seconds(exc: error.HTTPError) -> float | None:
        value = exc.headers.get("Retry-After") if exc.headers is not None else None
        if not value:
            return None
        try:
            return max(0.0, float(value))
        except ValueError:
            try:
                retry_at = parsedate_to_datetime(value)
                if retry_at.tzinfo is None:
                    retry_at = retry_at.replace(tzinfo=timezone.utc)
                now = datetime.now(retry_at.tzinfo)
                return max(0.0, (retry_at - now).total_seconds())
            except (TypeError, ValueError, OverflowError):
                return None

    def _get(self, endpoint: str, params: dict[str, Any]) -> bytes:
        params = {**params, "tool": "proteomics_pathway_report"}
        if self.email:
            params["email"] = self.email
        if self.api_key:
            params["api_key"] = self.api_key
        url = f"{PUBMED_EUTILS}/{endpoint}?{parse.urlencode(params)}"
        for attempt in range(len(NCBI_RETRY_DELAYS) + 1):
            self._throttle()
            try:
                with self._opener(url, timeout=self.timeout) as response:
                    return response.read()
            except error.HTTPError as exc:
                if (
                    exc.code not in NCBI_TRANSIENT_STATUS_CODES
                    or attempt == len(NCBI_RETRY_DELAYS)
                ):
                    suffix = (
                        " after 4 retries" if exc.code in NCBI_TRANSIENT_STATUS_CODES
                        else ""
                    )
                    raise RuntimeError(
                        f"PubMed request failed{suffix}: HTTP {exc.code} {exc.reason}"
                    ) from exc
                delay = float(NCBI_RETRY_DELAYS[attempt])
                retry_after = self._retry_after_seconds(exc)
                if retry_after is not None:
                    delay = max(delay, retry_after)
                retry_number = attempt + 1
                if exc.code == 429:
                    LOGGER.warning(
                        "PubMed rate limited (429); retry %d/%d in %g seconds.",
                        retry_number, len(NCBI_RETRY_DELAYS), delay,
                    )
                else:
                    LOGGER.warning(
                        "PubMed request failed with %d; retry %d/%d in %g seconds.",
                        exc.code, retry_number, len(NCBI_RETRY_DELAYS), delay,
                    )
                self._sleep(delay)
            except error.URLError as exc:
                raise RuntimeError(f"PubMed request failed: {exc}") from exc
        raise AssertionError("Unreachable PubMed retry state.")

    @staticmethod
    def build_query(theme: str, disease: str, today: date | None = None) -> str:
        today = today or date.today()
        start = today.replace(year=today.year - 5)
        safe_theme = theme.replace('"', "")
        safe_disease = disease.replace('"', "")
        return (
            f'("{safe_theme}"[Title/Abstract]) AND '
            f'("{safe_disease}"[Title/Abstract]) AND '
            f'("{start:%Y/%m/%d}"[Date - Publication] : "{today:%Y/%m/%d}"[Date - Publication])'
        )

    def search(self, theme: str, disease: str, *, max_results: int = 20,
               today: date | None = None) -> tuple[str, list[Study]]:
        query = self.build_query(theme, disease, today)
        root = ET.fromstring(self._get("esearch.fcgi", {
            "db": "pubmed", "term": query, "retmode": "xml", "retmax": max_results,
            "sort": "relevance",
        }))
        pmids = [node.text for node in root.findall(".//IdList/Id") if node.text]
        if not pmids:
            return query, []
        xml = self._get("efetch.fcgi", {
            "db": "pubmed", "id": ",".join(pmids), "retmode": "xml",
        })
        return query, self._parse_articles(xml)

    @staticmethod
    def _text(node: ET.Element | None) -> str:
        return "" if node is None else "".join(node.itertext()).strip()

    @classmethod
    def _parse_articles(cls, xml: bytes) -> list[Study]:
        studies: list[Study] = []
        for item in ET.fromstring(xml).findall(".//PubmedArticle"):
            medline = item.find("MedlineCitation")
            article = medline.find("Article") if medline is not None else None
            if medline is None or article is None:
                continue
            year_text = cls._text(article.find("Journal/JournalIssue/PubDate/Year"))
            if not year_text:
                year_text = cls._text(article.find("Journal/JournalIssue/PubDate/MedlineDate"))
            match = re.search(r"\b(19|20)\d{2}\b", year_text)
            abstract = " ".join(cls._text(x) for x in article.findall("Abstract/AbstractText"))
            authors = tuple(
                " ".join(filter(None, [cls._text(a.find("ForeName")), cls._text(a.find("LastName"))]))
                for a in article.findall("AuthorList/Author")
            )
            studies.append(Study(
                pmid=cls._text(medline.find("PMID")),
                title=cls._text(article.find("ArticleTitle")), abstract=abstract,
                journal=cls._text(article.find("Journal/Title")),
                year=int(match.group()) if match else None,
                authors=tuple(a for a in authors if a),
                publication_types=tuple(cls._text(x) for x in article.findall("PublicationTypeList/PublicationType")),
            ))
        return studies


class OpenAIResponsesClient:
    """Minimal Responses API adapter; avoids adding a runtime dependency."""

    def __init__(self, model: str, api_key: str | None = None, timeout: float = 120.0):
        self.model = model
        self.api_key = api_key or os.getenv("OPENAI_API_KEY")
        self.timeout = timeout
        if not self.api_key:
            raise ValueError("OPENAI_API_KEY is required to generate a report.")

    def json_response(self, instructions: str, payload: dict[str, Any]) -> dict[str, Any]:
        body = json.dumps({
            "model": self.model,
            "instructions": instructions,
            "input": json.dumps(payload, ensure_ascii=False),
            "text": {"format": {"type": "json_object"}},
        }).encode()
        req = request.Request(
            "https://api.openai.com/v1/responses", data=body, method="POST",
            headers={"Authorization": f"Bearer {self.api_key}", "Content-Type": "application/json"},
        )
        try:
            with request.urlopen(req, timeout=self.timeout) as response:
                result = json.loads(response.read())
        except error.HTTPError as exc:
            detail = exc.read().decode("utf-8", errors="replace")
            raise RuntimeError(f"OpenAI request failed ({exc.code}): {detail}") from exc
        text = result.get("output_text")
        if not text:
            chunks = [
                c.get("text", "") for output in result.get("output", [])
                for c in output.get("content", []) if c.get("type") == "output_text"
            ]
            text = "".join(chunks)
        if not text:
            raise RuntimeError("OpenAI response contained no output text.")
        return json.loads(text)


class GeminiClient:
    """Google Gemini adapter using the current official Google Gen AI SDK."""

    def __init__(
        self,
        model: str,
        api_key: str | None = None,
        *,
        client: Any | None = None,
        sleep: Callable[[float], None] = time.sleep,
    ):
        self.model = model
        self.api_key = api_key or os.getenv("GEMINI_API_KEY")
        if not self.api_key:
            raise ValueError("GEMINI_API_KEY is required to use the Gemini provider.")
        if client is None:
            try:
                from google import genai
            except ImportError as exc:
                raise RuntimeError(
                    "Gemini support requires the 'google-genai' package. "
                    "Install the project requirements first."
                ) from exc
            client = genai.Client(api_key=self.api_key)
        self._client = client
        self._sleep = sleep

    def generate(self, prompt: str) -> str:
        """Generate and return plain text from Gemini."""
        for attempt in range(len(GEMINI_RETRY_DELAYS) + 1):
            try:
                response = self._client.models.generate_content(
                    model=self.model,
                    contents=prompt,
                )
                break
            except Exception as exc:
                status_code = getattr(exc, "code", None)
                error_text = str(exc).lower()
                exhausted_quota = status_code == 429 and any(
                    marker in error_text
                    for marker in (
                        "per day", "daily quota", "project quota",
                        "quota has been exhausted", "quota exhausted",
                    )
                )
                if (
                    status_code not in GEMINI_TRANSIENT_STATUS_CODES
                    or exhausted_quota
                    or attempt == len(GEMINI_RETRY_DELAYS)
                ):
                    if exhausted_quota:
                        exc.add_note(
                            "Gemini reported daily/project quota exhaustion; "
                            "the request was not retried."
                        )
                    if (
                        status_code in GEMINI_TRANSIENT_STATUS_CODES
                        and attempt == len(GEMINI_RETRY_DELAYS)
                    ):
                        exc.add_note(
                            "Gemini request failed after 4 retries; "
                            f"last transient status was {status_code}."
                        )
                    raise
                delay = GEMINI_RETRY_DELAYS[attempt]
                retry_number = attempt + 1
                LOGGER.warning(
                    "Gemini request failed with %s; retry %d/%d in %d seconds.",
                    status_code,
                    retry_number,
                    len(GEMINI_RETRY_DELAYS),
                    delay,
                )
                self._sleep(delay)
        text = response.text
        if not text:
            raise RuntimeError("Gemini response contained no output text.")
        return str(text)

    def json_response(
        self, instructions: str, payload: dict[str, Any]
    ) -> dict[str, Any]:
        """Adapt Gemini plain-text generation to the agent's JSON interface."""
        prompt = "\n\n".join(
            [instructions, "Input JSON:", json.dumps(payload, ensure_ascii=False)]
        )
        text = self.generate(prompt).strip()
        if text.startswith("```"):
            lines = text.splitlines()
            if lines and lines[0].strip().lower() in {"```", "```json"}:
                lines = lines[1:]
            if lines and lines[-1].strip() == "```":
                lines = lines[:-1]
            text = "\n".join(lines).strip()
        try:
            result = json.loads(text)
        except json.JSONDecodeError as exc:
            raise RuntimeError("Gemini response was not valid JSON.") from exc
        if not isinstance(result, dict):
            raise RuntimeError("Gemini response JSON must be an object.")
        return result


class DatasetReportAgent:
    def __init__(
        self,
        llm: LLMClient,
        pubmed: PubMedClient,
        prompt: str,
        *,
        max_studies: int = 20,
        literature_batch_size: int = 4,
        abstract_max_chars: int = DEFAULT_ABSTRACT_MAX_CHARS,
        cache_dir: Path | None = None,
        cache_key: str | None = None,
        refresh_cache: bool = False,
    ):
        if literature_batch_size < 1:
            raise ValueError("literature_batch_size must be positive.")
        self.llm, self.pubmed, self.prompt, self.max_studies = llm, pubmed, prompt, max_studies
        self.literature_batch_size = literature_batch_size
        self.abstract_max_chars = abstract_max_chars
        self.cache_dir, self.cache_key = cache_dir, cache_key
        self.refresh_cache = refresh_cache
        self.request_counts = {
            "biological_themes": 0, "literature_assessment": 0,
            "report_writing": 0, "revision": 0,
        }
        self.last_run_context: dict[str, Any] | None = None

    def _llm_request(self, task: str, payload: dict[str, Any]) -> dict[str, Any]:
        count_key = {
            "name_biological_themes": "biological_themes",
            "assess_literature": "literature_assessment",
            "write_report": "report_writing",
            "revise_report": "revision",
        }[task]
        self.request_counts[count_key] += 1
        return self.llm.json_response(self.prompt, {"task": task, **payload})

    def request_summary(self, reviewer_calls: int = 0) -> str:
        total = sum(self.request_counts.values()) + reviewer_calls
        return "\n".join([
            "LLM request summary:",
            f"- biological theme calls: {self.request_counts['biological_themes']}",
            f"- literature assessment calls: {self.request_counts['literature_assessment']}",
            f"- report writing calls: {self.request_counts['report_writing']}",
            f"- reviewer calls: {reviewer_calls}",
            f"- revision calls: {self.request_counts['revision']}",
            f"- total LLM calls: {total}",
        ])

    @staticmethod
    def _render_report_response(
        report: dict[str, Any], clusters: list[dict[str, Any]], disease: str
    ) -> str:
        expected = {str(cluster["cluster_id"]) for cluster in clusters}
        sections = report.get("cluster_sections", [])
        returned = [str(item.get("cluster_id")) for item in sections]
        if len(returned) != len(set(returned)) or set(returned) != expected:
            raise ValueError("Report response must contain every cluster ID exactly once.")
        by_id = {
            str(item["cluster_id"]): str(item["markdown"]).strip()
            for item in sections
        }
        ordered = [by_id[str(cluster["cluster_id"])] for cluster in clusters]
        executive = str(report.get("executive_summary", "")).strip()
        final_interpretation = str(
            report.get("final_biological_interpretation", "")
        ).strip()
        return "\n\n".join([
            f"# Pathway Interpretation Report: {disease}",
            "## Executive Summary\n\n" + executive,
            *ordered,
            "## Final Biological Interpretation\n\n" + final_interpretation,
        ]).rstrip() + "\n"

    @staticmethod
    def _input_hash(value: Any) -> str:
        serialized = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
        return hashlib.sha256(serialized.encode("utf-8")).hexdigest()

    def _cache_path(self, kind: str) -> Path | None:
        if self.cache_dir is None or self.cache_key is None:
            return None
        return self.cache_dir / f"{kind}_{self.cache_key}.json"

    def _read_cache(self, kind: str, inputs: Any) -> Any | None:
        path = self._cache_path(kind)
        if path is None or self.refresh_cache or not path.exists():
            return None
        try:
            cached = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            return None
        if not isinstance(cached, dict):
            return None
        if cached.get("input_hash") != self._input_hash(inputs):
            return None
        return cached.get("result")

    def _write_cache(self, kind: str, inputs: Any, result: Any) -> None:
        path = self._cache_path(kind)
        if path is None:
            return
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_suffix(path.suffix + ".tmp")
        temporary.write_text(json.dumps({
            "input_hash": self._input_hash(inputs), "result": result,
        }, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        temporary.replace(path)

    @staticmethod
    def _theme_input(cluster: dict[str, Any]) -> dict[str, Any]:
        return {
            "cluster_id": cluster["cluster_id"],
            "representative_pathway": cluster["representative_pathway"],
            "member_pathways": cluster["member_pathways"],
            "supporting_genes": cluster["genes"],
        }

    def build_search_themes(
        self, evidence: dict[str, Any], disease: str
    ) -> dict[str, str]:
        """Generate and validate the search theme for every retained cluster."""
        clusters = evidence["clusters"]
        inputs = {
            "prompt": self.prompt, "disease": disease,
            "clusters": [self._theme_input(c) for c in clusters],
        }
        themes_result = self._read_cache("biological_themes", inputs)
        if themes_result is None:
            themes_result = self._llm_request("name_biological_themes", {
                "disease": disease, "clusters": inputs["clusters"],
            })
        theme_items = themes_result.get("themes", [])
        theme_ids = [str(x.get("cluster_id")) for x in theme_items]
        if len(theme_ids) != len(set(theme_ids)):
            raise ValueError("Theme response duplicated a cluster ID.")
        themes = {str(x["cluster_id"]): str(x["biological_theme"]) for x in theme_items}
        expected = {str(c["cluster_id"]) for c in clusters}
        if set(themes) != expected:
            raise ValueError("Theme response must contain every cluster ID exactly once.")
        self._write_cache("biological_themes", inputs, themes_result)
        return themes

    def _study_payload(self, study: Study) -> dict[str, Any]:
        value = asdict(study) | {
            "is_review": study.is_review,
            "is_non_primary_publication": study.is_non_primary_publication,
        }
        value["abstract"] = study.abstract[:self.abstract_max_chars]
        return value

    @staticmethod
    def _validate_assessment(
        cluster_id: str, assessment: dict[str, Any], studies: list[Study]
    ) -> None:
        annotations = assessment.get("study_annotations", [])
        retrieved_pmids = {study.pmid for study in studies}
        annotated_pmids = [str(item.get("pmid")) for item in annotations]
        if len(annotated_pmids) != len(set(annotated_pmids)):
            raise ValueError(f"Literature response duplicated a PMID for cluster {cluster_id}.")
        if set(annotated_pmids) != retrieved_pmids:
            raise ValueError(
                "Literature response must assess every retrieved PMID exactly once "
                f"for cluster {cluster_id}."
            )
        studies_by_pmid = {study.pmid: study for study in studies}
        invalid_primary = [
            str(item["pmid"]) for item in annotations
            if item.get("is_primary_study")
            and studies_by_pmid[str(item["pmid"])].is_non_primary_publication
        ]
        if invalid_primary:
            raise ValueError(
                "Literature response classified non-primary publication types "
                f"as primary for cluster {cluster_id}: {invalid_primary}"
            )

    def _assess_literature_batches(
        self, disease: str, retrievals: list[dict[str, Any]]
    ) -> dict[str, dict[str, Any]]:
        model_clusters = [
            {key: value for key, value in item.items() if key != "_studies"}
            for item in retrievals
        ]
        inputs = {"prompt": self.prompt, "disease": disease, "clusters": model_clusters}
        cached = self._read_cache("literature_assessment", inputs)
        if cached is not None:
            assessment_items = cached.get("cluster_assessments", [])
        else:
            assessment_items = []
            for start in range(0, len(model_clusters), self.literature_batch_size):
                batch = model_clusters[start:start + self.literature_batch_size]
                response = self._llm_request("assess_literature", {
                    "disease": disease, "clusters": batch,
                })
                assessment_items.extend(response.get("cluster_assessments", []))
        ids = [str(item.get("cluster_id")) for item in assessment_items]
        expected = {str(item["cluster_id"]) for item in retrievals}
        if len(ids) != len(set(ids)) or set(ids) != expected:
            raise ValueError("Literature response must contain every cluster ID exactly once.")
        assessments = {str(item["cluster_id"]): item for item in assessment_items}
        for retrieval in retrievals:
            cluster_id = str(retrieval["cluster_id"])
            self._validate_assessment(cluster_id, assessments[cluster_id], retrieval["_studies"])
        if cached is None:
            self._write_cache(
                "literature_assessment", inputs,
                {"cluster_assessments": assessment_items},
            )
        return assessments

    def preview_queries(
        self, evidence: dict[str, Any], disease: str
    ) -> list[dict[str, str]]:
        """Build search queries without calling PubMed or generating a report."""
        self.request_counts = {key: 0 for key in self.request_counts}
        themes = self.build_search_themes(evidence, disease)
        return [
            {
                "cluster_id": str(cluster["cluster_id"]),
                "representative_pathway": str(cluster["representative_pathway"]["name"]),
                "biological_search_theme": themes[str(cluster["cluster_id"])],
                "pubmed_query": self.pubmed.build_query(
                    themes[str(cluster["cluster_id"])], disease
                ),
            }
            for cluster in evidence["clusters"]
        ]

    def generate(self, evidence: dict[str, Any], disease: str) -> str:
        self.request_counts = {key: 0 for key in self.request_counts}
        clusters = evidence["clusters"]
        themes = self.build_search_themes(evidence, disease)
        retrievals: list[dict[str, Any]] = []
        for cluster in clusters:
            cluster_id = str(cluster["cluster_id"])
            query, studies = self.pubmed.search(themes[cluster_id], disease, max_results=self.max_studies)
            retrievals.append({
                "cluster_id": cluster_id,
                "biological_theme": themes[cluster_id],
                "pubmed_query": query,
                "pubmed_records": [self._study_payload(study) for study in studies],
                "_studies": studies,
            })
        assessments = self._assess_literature_batches(disease, retrievals)

        packets: list[dict[str, Any]] = []
        for cluster, retrieval in zip(clusters, retrievals):
            cluster_id = str(cluster["cluster_id"])
            literature = assessments[cluster_id]
            annotations = literature.get("study_annotations", [])
            primary_count = sum(
                bool(a.get("is_primary_study") and a.get("direct_disease_support"))
                for a in annotations
            )
            counts = {int(k): int(v) for k, v in cluster["level_counts"].items()}
            packets.append({
                "cluster": cluster, "biological_theme": themes[cluster_id],
                "dataset_evidence_strength": dataset_evidence_strength(cluster),
                "total_proteins": sum(counts.values()),
                "levels_1_3": sum(counts.get(x, 0) for x in (1, 2, 3)),
                "levels_4_5": sum(counts.get(x, 0) for x in (4, 5)),
                "pubmed_query": retrieval["pubmed_query"],
                "retrieved_studies": [asdict(s) for s in retrieval["_studies"]],
                "literature_assessment": literature,
                "relevant_primary_study_count": primary_count,
                "disease_relevance": disease_relevance(annotations),
            })

        report = self._llm_request("write_report", {
            "disease": disease,
            "comparison": evidence.get("metadata", {}).get("comparison", "unknown"),
            "clusters": packets,
        })
        self.last_run_context = {
            "disease": disease,
            "pathway_evidence": evidence,
            "biological_themes": themes,
            "retrieved_pubmed": [
                {key: value for key, value in retrieval.items() if key != "_studies"}
                for retrieval in retrievals
            ],
            "literature_assessments": assessments,
            "deterministic_clusters": packets,
        }
        return self._render_report_response(report, clusters, disease)

    def revise_report(
        self, draft_report: str, reviewer_errors: list[dict[str, Any]]
    ) -> str:
        """Ask the report writer to fix only reviewer-identified errors."""
        if self.last_run_context is None:
            raise RuntimeError("A report must be generated before it can be revised.")
        response = self._llm_request("revise_report", {
            "disease": self.last_run_context["disease"],
            "clusters": self.last_run_context["deterministic_clusters"],
            "original_draft_report": draft_report,
            "reviewer_errors": reviewer_errors,
        })
        return self._render_report_response(
            response,
            self.last_run_context["pathway_evidence"]["clusters"],
            self.last_run_context["disease"],
        )
