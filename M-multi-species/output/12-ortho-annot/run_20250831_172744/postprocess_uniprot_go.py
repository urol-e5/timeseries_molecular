import sys, os, time, json, math
import pandas as pd
import requests
from pathlib import Path

BLAST_TSV = sys.argv[1]
OUTDIR = sys.argv[2]

os.makedirs(OUTDIR, exist_ok=True)
hits = pd.read_csv(BLAST_TSV, sep='\t', header=None,
                   names=['query','accession','pident','length','evalue','bitscore','title'])

# Normalize accession (handle formats like "sp|P69907|HBA_PANTR") for UniProt REST API
def _norm_acc(a: str):
    if not isinstance(a, str):
        return a
    if '|' in a:
        parts = a.split('|')
        if len(parts) >= 3 and parts[0] in {'sp','tr'}:
            return parts[1]
    return a
hits['accession'] = hits['accession'].apply(_norm_acc)
hits['query'] = hits['query'].astype(str)

# Deduplicate accessions for API efficiency
accs = sorted(hits['accession'].dropna().unique().tolist())
# UniProt REST batch: use query with OR list; fields: accession, id, protein_name, organism_name, go_id, go_p, go_c, go_f
# See UniProt REST docs. We will page in chunks to be safe.
# https://www.uniprot.org/help/api_queries ; https://www.uniprot.org/help/return_fields
def fetch_uniprot_batch(acc_batch):
    """Fetch UniProt annotations for a list of accessions, splitting recursively if query causes 400."""
    if not acc_batch:
        return pd.DataFrame()
    # Base case: single accession -> always try directly
    if len(acc_batch) == 1:
        query = f"accession:{acc_batch[0]}"
    else:
        query = " OR ".join([f"accession:{a}" for a in acc_batch])
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": query,
        "format": "tsv",
        "fields": "accession,id,reviewed,protein_name,organism_name,go_id,go_p,go_c,go_f"
    }
    try:
        r = requests.get(url, params=params, timeout=60)
        r.raise_for_status()
    except requests.HTTPError as e:
        # If too large / bad request and we have multiple accs, split and recurse
        if r.status_code == 400 and len(acc_batch) > 1:
            mid = len(acc_batch)//2
            left = fetch_uniprot_batch(acc_batch[:mid])
            right = fetch_uniprot_batch(acc_batch[mid:])
            return pd.concat([left, right], ignore_index=True)
        raise
    if not r.text.strip():
        return pd.DataFrame(columns=["accession","id","reviewed","protein_name","organism","go_ids","go_bp","go_cc","go_mf"])
    from io import StringIO
    df = pd.read_csv(StringIO(r.text), sep='\t')
    # normalize column names
    rename_map = {
        "Entry":"accession",
        "Entry Name":"id",
        "Status":"reviewed",
        "Protein names":"protein_name",
        "Organism":"organism",
        "Gene Ontology IDs":"go_ids",
        "Gene Ontology (biological process)":"go_bp",
        "Gene Ontology (cellular component)":"go_cc",
        "Gene Ontology (molecular function)":"go_mf"
    }
    df = df.rename(columns={k:v for k,v in rename_map.items() if k in df.columns})
    return df

frames = []
CHUNK = 50  # smaller initial chunk to avoid long URLs
for i in range(0, len(accs), CHUNK):
    batch = accs[i:i+CHUNK]
    frames.append(fetch_uniprot_batch(batch))
    time.sleep(0.1)

if frames:
    u = pd.concat(frames, ignore_index=True)
else:
    # Ensure merge does not fail if there are no hits
    u = pd.DataFrame(columns=['accession','id','reviewed','protein_name','organism','go_ids','go_bp','go_cc','go_mf'])
# Join back to BLAST best hits
merged = hits.merge(u, how='left', on='accession')

# Write intermediate with full GO
full_path = Path(OUTDIR) / "annotation_full_go.tsv"
merged.to_csv(full_path, sep='\t', index=False)

# ---- Map full GO -> GO-Slim (generic) using GOATOOLS ----
# inputs: go-basic.obo + goslim_generic.obo
from goatools.obo_parser import GODag
from collections import defaultdict

obo_basic = Path(OUTDIR) / "go-basic.obo"
obo_slim  = Path(OUTDIR) / "goslim_generic.obo"

# download if missing
import urllib.request
if not obo_basic.exists():
    urllib.request.urlretrieve("http://purl.obolibrary.org/obo/go/go-basic.obo", obo_basic.as_posix())
if not obo_slim.exists():
    urllib.request.urlretrieve("https://go.princeton.edu/GOTermMapper/goSlimFiles/goslim_generic.obo", obo_slim.as_posix())

obodag = GODag(obo_basic.as_posix(), optional_attrs={'relationship'})
slimdag = GODag(obo_slim.as_posix(), optional_attrs={'relationship'})

# Build query: list of GO terms per gene (merge bp/cc/mf) as GO IDs
def split_goids(s):
    if isinstance(s, float) and math.isnan(s): return []
    # UniProt 'go_ids' is semi-colon separated IDs (e.g., "GO:0005524; GO:0004672")
    return [x.strip() for x in str(s).split(';') if x.strip().startswith('GO:')]

gene2gos = defaultdict(set)
for _, row in merged.iterrows():
    gs = set(split_goids(row.get('go_ids', '')))
    if gs:
        gene2gos[row['query']].update(gs)

# Map to slim: returns counts per slim term; we also want per-gene mapping
# We'll perform per-gene slim mapping using GOATOOLS internal search
def _ancestors(goid, dag):
    term = dag.get(goid)
    if term is None:
        return []
    seen = set()
    stack = [term]
    while stack:
        t = stack.pop()
        if t.id in seen:
            continue
        seen.add(t.id)
        for p in t.parents:
            stack.append(p)
    return seen

slim_set = {t.id for t in slimdag.values()}
slim_map_rows = []
for gene, gos in gene2gos.items():
    slim_hits = set()
    for go in gos:
        for anc in _ancestors(go, obodag):
            if anc in slim_set:
                slim_hits.add(anc)
    slim_map_rows.append({
        "query": gene,
        "goslim_ids": "; ".join(sorted(slim_hits))
    })
if slim_map_rows:
    slim_df = pd.DataFrame(slim_map_rows)
    # Add names for slim terms
    def name_of(goid):
        n = slimdag.get(goid)
        return n.name if n is not None else ""
    slim_df['goslim_names'] = slim_df['goslim_ids'].apply(
        lambda s: "; ".join([name_of(x) for x in s.split('; ')]) if s else ""
    )
else:
    # Ensure required columns exist for merge
    slim_df = pd.DataFrame(columns=["query","goslim_ids","goslim_names"])
    slim_df = slim_df.astype({"query":"string","goslim_ids":"string","goslim_names":"string"})

final = merged.merge(slim_df, on='query', how='left')
final['query'] = final['query'].astype(str)

# Tidy final columns
desired_cols = [
    "query","accession","id","reviewed","protein_name","organism",
    "pident","length","evalue","bitscore","title",
    "go_ids","go_bp","go_cc","go_mf",
    "goslim_ids","goslim_names"
]
for c in desired_cols:
    if c not in final.columns:
        final[c] = pd.NA
final = final[desired_cols]

final.to_csv(Path(OUTDIR) / "annotation_with_goslim.tsv", sep='\t', index=False)
print(f"[DONE] Wrote:\n - {full_path}\n - {Path(OUTDIR) / 'annotation_with_goslim.tsv'}")

# Generate summary report
blast_hits = len(final[final['accession'].notna()])
go_matches = len(final[final['go_ids'].notna() & (final['go_ids'] != '')])
goslim_matches = len(final[final['goslim_ids'].notna() & (final['goslim_ids'] != '')])
total_seqs = len(final)

# GO-Slim category counts for visualization (Biological Process only)
goslim_counts = {}
if goslim_matches > 0:
    for _, row in final[final['goslim_names'].notna() & (final['goslim_names'] != '')].iterrows():
        terms = [x.strip() for x in str(row['goslim_names']).split(';') if x.strip()]
        for term in terms:
            goslim_counts[term] = goslim_counts.get(term, 0) + 1

# Filter GO-Slim counts for Biological Process terms only
bp_goslim_counts = {}
if goslim_counts:
    # Get corresponding GO IDs for filtering by namespace
    for _, row in final[final['goslim_ids'].notna() & (final['goslim_ids'] != '') & 
                       final['goslim_names'].notna() & (final['goslim_names'] != '')].iterrows():
        slim_ids = [x.strip() for x in str(row['goslim_ids']).split(';') if x.strip()]
        slim_names = [x.strip() for x in str(row['goslim_names']).split(';') if x.strip()]
        
        for slim_id, slim_name in zip(slim_ids, slim_names):
            if slim_id in slimdag:
                term_obj = slimdag[slim_id]
                # Filter for biological process namespace
                if hasattr(term_obj, 'namespace') and term_obj.namespace == 'biological_process':
                    bp_goslim_counts[slim_name] = bp_goslim_counts.get(slim_name, 0) + 1

# Use BP counts for visualization, fall back to all counts if no BP terms
chart_counts = bp_goslim_counts if bp_goslim_counts else goslim_counts

# Write summary stats for main script
summary_stats = {
    'total_sequences': total_seqs,
    'blast_hits': blast_hits,
    'go_matches': go_matches,
    'goslim_matches': goslim_matches,
    'goslim_counts': chart_counts  # Use filtered BP counts for chart
}

import json
with open(Path(OUTDIR) / "summary_stats.json", 'w') as f:
    json.dump(summary_stats, f)

print(f"[INFO] Summary stats: {blast_hits}/{total_seqs} BLAST hits, {go_matches} GO annotations, {goslim_matches} GO-Slim mappings")
