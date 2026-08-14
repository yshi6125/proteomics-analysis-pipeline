"""CLI entry point for disease-specific pathway interpretation reports."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import sys

if __package__ is None:  # Support ``python src/generate_pathway_report.py``.
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.agents.dataset_report_agent import (
    DatasetReportAgent, GeminiClient, LLMClient, OpenAIResponsesClient,
    PubMedClient, load_pathway_evidence,
)
from src.agents.scientific_reviewer import ScientificReviewerAgent


PROJECT_ROOT = Path(__file__).resolve().parent.parent


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--disease", required=True, help="Disease for the literature review")
    parser.add_argument("--comparison", required=True, help="Comparison suffix, e.g. DM_vs_NDM")
    parser.add_argument("--branch", choices=("primary", "exploratory", "both"), default="both")
    parser.add_argument(
        "--provider", choices=("openai", "gemini"), default="openai",
        help="LLM provider (default: openai)",
    )
    parser.add_argument("--model", default="gpt-5-mini", help="OpenAI model name")
    parser.add_argument(
        "--reviewer-model",
        help="Reviewer model; defaults to --model using the same provider",
    )
    parser.add_argument(
        "--skip-review", action="store_true",
        help="Skip scientific review for development/debugging",
    )
    parser.add_argument("--max-studies", type=int, default=20)
    parser.add_argument(
        "--literature-batch-size", type=int, default=4,
        help="Clusters per literature-assessment LLM request (default: 4)",
    )
    parser.add_argument(
        "--refresh-cache", action="store_true",
        help="Ignore and replace cached themes and literature assessments",
    )
    parser.add_argument(
        "--preview-queries", action="store_true",
        help="Print biological themes and PubMed queries without retrieval or report generation",
    )
    parser.add_argument("--evidence", type=Path, help="Explicit evidence JSON (one branch only)")
    parser.add_argument("--output", type=Path, help="Explicit output Markdown (one branch only)")
    return parser.parse_args()


def build_llm_client(provider: str, model: str) -> LLMClient:
    """Construct the selected provider client."""
    if provider == "openai":
        return OpenAIResponsesClient(model)
    if provider == "gemini":
        return GeminiClient(model)
    raise ValueError(f"Unsupported LLM provider: {provider}")


def cache_key(comparison: str, disease: str) -> str:
    disease_slug = re.sub(r"[^a-z0-9]+", "_", disease.lower()).strip("_")
    return f"{comparison}_{disease_slug or 'disease'}"


def main() -> None:
    args = parse_args()
    if (args.evidence or args.output) and args.branch == "both":
        raise SystemExit("--evidence/--output require --branch primary or exploratory")
    if args.max_studies < 1:
        raise SystemExit("--max-studies must be positive")
    if args.literature_batch_size < 1:
        raise SystemExit("--literature-batch-size must be positive")
    prompt = (PROJECT_ROOT / "prompts" / "dataset_report_prompt.txt").read_text(encoding="utf-8")
    reviewer_prompt = (
        PROJECT_ROOT / "prompts" / "scientific_reviewer_prompt.txt"
    ).read_text(encoding="utf-8")
    llm = build_llm_client(args.provider, args.model)
    reviewer_model = args.reviewer_model or args.model
    reviewer_llm = (
        llm if reviewer_model == args.model
        else build_llm_client(args.provider, reviewer_model)
    )
    branches = ("primary", "exploratory") if args.branch == "both" else (args.branch,)
    for branch in branches:
        report_dir = PROJECT_ROOT / "results" / "pathway_reports" / branch
        agent = DatasetReportAgent(
            llm, PubMedClient(), prompt,
            max_studies=args.max_studies,
            literature_batch_size=args.literature_batch_size,
            cache_dir=report_dir / "cache",
            cache_key=cache_key(args.comparison, args.disease),
            refresh_cache=args.refresh_cache,
        )
        evidence_path = args.evidence or (
            PROJECT_ROOT / "results" / "pathway_consolidation" / branch
            / f"pathway_evidence_{args.comparison}.json"
        )
        evidence = load_pathway_evidence(evidence_path)
        if args.preview_queries:
            for item in agent.preview_queries(evidence, args.disease):
                print(f"cluster_id: {item['cluster_id']}")
                print(f"representative pathway: {item['representative_pathway']}")
                print(f"biological search theme: {item['biological_search_theme']}")
                print(f"PubMed query: {item['pubmed_query']}")
                print()
            print(agent.request_summary())
            continue
        output_path = args.output or (
            report_dir / f"dataset_report_{args.comparison}.md"
        )
        report = agent.generate(evidence, args.disease)
        reviewer_calls = 0
        if not args.skip_review:
            if agent.last_run_context is None:
                raise RuntimeError("Report generation did not retain reviewer context.")
            reviewer = ScientificReviewerAgent(
                reviewer_llm, reviewer_prompt,
                provider=args.provider, model=reviewer_model,
                cache_dir=report_dir / "cache",
                cache_key=cache_key(args.comparison, args.disease),
                refresh_cache=args.refresh_cache,
            )
            outcome = reviewer.run(
                report, agent.last_run_context, agent.revise_report
            )
            report = outcome.final_report
            reviewer_calls = reviewer.request_count
            review_dir = report_dir / "review"
            review_dir.mkdir(parents=True, exist_ok=True)
            metadata = {
                "reviewer_enabled": True,
                "review_status": outcome.review_status,
                "review_rounds_used": outcome.review_rounds_used,
                "reviewer_model": reviewer_model,
                "provider": args.provider,
                "final_report_model": args.model,
                "generated_at": datetime.now(timezone.utc).isoformat(),
            }
            review_payload = {
                "metadata": metadata,
                "status": outcome.final_review["status"],
                "errors": outcome.final_review["errors"],
            }
            review_path = review_dir / f"review_result_{args.comparison}.json"
            review_path.write_text(
                json.dumps(review_payload, indent=2, ensure_ascii=False) + "\n",
                encoding="utf-8",
            )
            unresolved_path = (
                review_dir / f"unresolved_review_issues_{args.comparison}.json"
            )
            if outcome.review_status == "manual_review_required":
                unresolved_path.write_text(
                    json.dumps(review_payload, indent=2, ensure_ascii=False) + "\n",
                    encoding="utf-8",
                )
            else:
                unresolved_path.unlink(missing_ok=True)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(report, encoding="utf-8")
        print(f"Saved {branch} report: {output_path}")
        print(agent.request_summary(reviewer_calls))


if __name__ == "__main__":
    main()
