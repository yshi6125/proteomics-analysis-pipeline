"""Historical/reference reproduction only. Production commands use pathway_protein and pathway_comparison."""
import argparse
from pathlib import Path

from src.pathway_analysis import run_pathway_analysis
from src.pathway_consolidation import consolidate_pathways
from src.pathway_measured_ora import run_measured_ora
from src.pathway_comparison import run_comparison


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage',choices=['legacy-analysis','measured-analysis','comparison'],required=True)
    parser.add_argument('--comparison',default='DM_vs_NDM')
    parser.add_argument('--gene-sets',type=Path)
    parser.add_argument('--a-mode',choices=['legacy','measured'],default='legacy',help='Reference comparison only')
    parser.add_argument('--a-branch',choices=['primary','exploratory'],default='primary',help='Reference comparison only')
    args=parser.parse_args()
    print('REFERENCE ONLY: historical alternatives, outside the production workflow.')
    if args.stage=='legacy-analysis':
        run_pathway_analysis(comparison=args.comparison)
        consolidate_pathways(comparison=args.comparison)
    elif args.stage=='measured-analysis':
        if args.gene_sets is None:
            parser.error('--gene-sets is required for historical measured-analysis')
        run_measured_ora(args.comparison,gene_sets_path=args.gene_sets)
    else:
        run_comparison(args.comparison,a_mode=args.a_mode,branch=args.a_branch,gene_sets_path=args.gene_sets)


if __name__=='__main__':
    main()
