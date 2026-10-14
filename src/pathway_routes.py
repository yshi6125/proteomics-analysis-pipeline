"""Run production Route A, Route B, or both independently. Production has a single Route A scope."""
import argparse
import json
from pathlib import Path

from src.pathway_protein import run_protein_centric, default_snapshot
from src.pathway_centric import run_pathway_centric, DEFAULT_LIBRARIES


def run_routes(route='protein', comparison='DM_vs_NDM', *, gene_sets_path=None, **pathway_options):
    if route not in {'protein','pathway','both'}:
        raise ValueError('route must be protein, pathway, or both')
    outputs={}
    snapshot=Path(gene_sets_path) if gene_sets_path is not None else None
    if route in {'protein','both'}:
        snapshot=snapshot or default_snapshot(comparison)
        outputs['protein']=run_protein_centric(comparison=comparison,gene_sets_path=snapshot)
    if route in {'pathway','both'}:
        options=dict(pathway_options)
        if snapshot is not None:
            if 'gene_sets' in options:
                raise ValueError('Supply gene_sets_path or gene_sets, not both')
            options['gene_sets']=json.loads(snapshot.read_text())
        outputs['pathway']=run_pathway_centric(comparison=comparison,**options)
    return outputs


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--route',choices=['protein','pathway','both'],default='protein',
                        help='protein = Route A; pathway = Route B')
    parser.add_argument('--comparison',default='DM_vs_NDM')
    parser.add_argument('--gene-sets',type=Path,help='Pinned local JSON; both routes use this same snapshot')
    parser.add_argument('--libraries',nargs='+',default=list(DEFAULT_LIBRARIES),help='Route B live-library option only')
    parser.add_argument('--organism',default='Human',help='Route B option only')
    parser.add_argument('--min-size',type=int,default=5,help='Route B option only')
    parser.add_argument('--max-size',type=int,default=500,help='Route B option only')
    parser.add_argument('--permutations',type=int,default=100000,help='Route B option only (default: 100000)')
    parser.add_argument('--seed',type=int,default=42,help='Route B option only')
    parser.add_argument('--q-cutoff',type=float,default=.05,help='Route B option only')
    parser.add_argument('--p-cutoff',type=float,default=1.,help='Route B option only')
    parser.add_argument('--jaccard-threshold',type=float,default=.25,help='Route B option only; use pathway_protein for Route A parameters')
    parser.add_argument('--louvain-seed',type=int,default=42,help='Route B option only')
    args=vars(parser.parse_args())
    args['gene_sets_path']=args.pop('gene_sets')
    outputs=run_routes(**args)
    if 'pathway' in outputs:
        print(json.dumps(outputs['pathway'][1],indent=2))


if __name__=='__main__':
    main()
