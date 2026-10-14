from pathlib import Path
import re,json,csv,collections
root=Path.cwd()
pattern=re.compile(r'measured_exploratory|exploratory_all_candidates|legacy|exploratory|primary',re.I)
excluded={'docs/CLEANUP_SEARCH_AUDIT.tsv','docs/CLEANUP_SEARCH_SUMMARY.json','docs/CLEANUP_REPORT.md'}
rows=[]
production_modules={'src/pathway_protein.py','src/pathway_ora.py','src/protein_evidence.py','src/pathway_redundancy.py','src/pathway_comparison.py','src/pathway_centric.py','src/pathway_routes.py','src/pathway_review_io.py'}
def category(rel,line):
    if rel in production_modules:return 'ACTIVE','Production module match requiring review'
    if rel.startswith('archive/') or 'archive/route_a_legacy' in line:return 'HISTORICAL/ARCHIVED','Archived record or explicit archive link/path; no production workflow choice'
    if rel.startswith('tests/') or rel.startswith('docs/') or rel.startswith('scripts/'):return 'TEST/REFERENCE','Validation or explanatory documentation, not an active alternative'
    if rel=='README.md':return 'TEST/REFERENCE','Archive/cleanup explanation, not a production option'
    if rel.startswith(('results/pathway_protein/','results/pathway_comparison/')) and 'historical_method_reference' in line:
        return 'HISTORICAL/ARCHIVED','Immutable recorded generation-time method identity; archive manifest resolves original path'
    if 'differential_abundance' in rel:return 'ACTIVE','Shared protein statistics, annotations or tier diagnostics; not a retired pathway workflow'
    if rel.startswith('src/agents/'):
        return 'ACTIVE','Unchanged literature primary-study or dataset-strength classification, not a Route A branch'
    if rel.startswith(('data/', 'results/')):
        return 'ACTIVE','Biological protein/pathway annotation in shared data or production output; not a retired workflow choice'
    if rel=='src/normalization_diagnostics.py':
        return 'ACTIVE','Shared diagnostic interpretation; not a Route A workflow choice'
    if rel=='.gitignore':return 'TEST/REFERENCE','Ignore policy for archived/generated artifacts'
    return 'TEST/REFERENCE','Nonproduction reference record; inspected occurrence'
for p in sorted(root.rglob('*')):
    if not p.is_file():continue
    rel=p.relative_to(root).as_posix()
    if rel in excluded or any(part in {'.git','__pycache__','.pytest_cache','.mypy_cache','.venv'} for part in p.relative_to(root).parts):continue
    if p.suffix.lower() in {'.png','.jpg','.jpeg','.pdf','.xlsx','.pkl','.pyc','.ds_store'} or p.name=='.DS_Store':continue
    try:text=p.read_text()
    except (UnicodeError,OSError):continue
    for match in pattern.finditer(rel):
        label,why=category(rel,rel)
        rows.append([rel,0,match.start()+1,match.group(),label,why,'FILE PATH'])
    for number,line in enumerate(text.splitlines(),1):
        for match in pattern.finditer(line):
            label,why=category(rel,line)
            context=line[max(0,match.start()-65):match.end()+85].replace('\t',' ')
            rows.append([rel,number,match.start()+1,match.group(),label,why,context])
with (root/'docs/CLEANUP_SEARCH_AUDIT.tsv').open('w') as out:
    writer=csv.writer(out,delimiter='\t');writer.writerow(['path','line','column','match','classification','rationale','context']);writer.writerows(rows)
summary=dict(occurrences=len(rows),classifications=dict(collections.Counter(r[4] for r in rows)),
             active_files=sorted({r[0] for r in rows if r[4]=='ACTIVE'}),
             terms=['legacy','exploratory','primary','measured_exploratory','exploratory_all_candidates'],
             scope='All readable text and paths, including archive/dotfiles; excludes Git internals, caches, binaries and generated cleanup audit/report',
             exclusions=sorted(excluded),
             unintended_active_route_a_workflow_occurrences=sum(r[0] in production_modules for r in rows))
(root/'docs/CLEANUP_SEARCH_SUMMARY.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
print(json.dumps(summary,indent=2))
