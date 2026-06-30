#!/usr/bin/env python3
"""
Master summary table of the 7 species-novel MAGs — the data backbone for the
"novel bacteria" supplementary section that applies the Wang/Gao/Wells supinfo
(emi15486-sup-0001-supinfo.docx) to this project's novel organisms.

Merges, per novel MAG:
  • taxonomy + species-level novelty  (Table S1 xlsx + gtdbtk.bac120.summary.tsv)
  • genome quality / assembly stats   (Table S1 xlsx: CheckM2 compl/contam, size, N50, GC)
  • relative-abundance trajectory      (Table S1 xlsx CAN_1..5; + trend)
  • N-cycle role + denitrification/DNRA gene inventory + nosZ clade/subclade
                                       (gene_ab_figure/data/*)
  • carbon/storage/polyP/pili/B12/biotin functional-marker completeness
                                       (data/novel_functional_ko_long.tsv, this pipeline)

Outputs (./data): novel_bacteria_summary.tsv, novel_bacteria_summary.md
"""
import csv, os
import openpyxl

WORK = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(WORK)
DATA = os.path.join(WORK, 'data')

NOVEL = ['coasm_bin.55', 'CAN_2_bin.6', 'CAN_1_bin.210', 'CAN_1_bin.77',
         'CAN_5_bin.112', 'coasm_bin.260', 'coasm_bin.481']
SAMPLES = ['CAN_1', 'CAN_2', 'CAN_3', 'CAN_4', 'CAN_5']

def load_tsv(p):
    with open(p) as f:
        return list(csv.DictReader(f, delimiter='\t'))

# ── Table S1 (xlsx): quality / taxonomy / abundance ──────────────────────────
wb = openpyxl.load_workbook(os.path.join(ROOT, 'MAG_target_gene_filled.xlsx'), data_only=True)
ws = wb['Table S1_MAG quality']
# columns: 0 bin,1 iterID,2 phylum,3 class,4 order,5 family,6 genus,7 compl,
#          8 contam,9 size_bp,10 contigs,11 GC,12 N50,13-17 abund CAN1..5,18 lineage
S1 = {}
for r in ws.iter_rows(values_only=True):
    if r[0] in NOVEL:
        S1[r[0]] = r

# ── GTDB-Tk: RED + classification method ─────────────────────────────────────
red, method = {}, {}
for r in load_tsv(os.path.join(ROOT, 'gtdbtk.bac120.summary.tsv')):
    if r['user_genome'] in NOVEL:
        red[r['user_genome']] = r['red_value']
        method[r['user_genome']] = r['classification_method']

# ── denitrification inventory + nosZ clade/subclade + N-cycle role ───────────
cp = {r['MAG']: r for r in load_tsv(os.path.join(DATA, 'gene_copy_per_mag.tsv'))}
role = {r['MAG']: r['classification'] for r in load_tsv(os.path.join(DATA, 'selected_mags.tsv'))}

# denitrification step -> member columns (operon copy = max over subunits)
DENIT = [('nar', ['narG', 'narH', 'narI']), ('nap', ['napA', 'napB']),
         ('nirS', ['nirS']), ('nirK', ['nirK']), ('nirB/D(DNRA)', ['nirB', 'nirD']),
         ('nor', ['norB', 'norC']), ('nosZ', ['nosZ_cladeI', 'nosZ_cladeII'])]

def denit_str(m):
    r = cp[m]
    parts = []
    for label, members in DENIT:
        n = max(int(float(r.get(c, 0))) for c in members)
        if n > 0:
            parts.append(f'{label}' + (f'×{n}' if n > 1 else ''))
    return ', '.join(parts) if parts else '—'

def nosz_str(m):
    r = cp[m]
    cI, cII = int(float(r['nosZ_cladeI'])), int(float(r['nosZ_cladeII']))
    sub = r.get('nosZ_II_subclade', '').strip()
    if cI:
        return 'Clade I'
    if cII:
        return f'Clade II ({sub})' if sub else 'Clade II'
    return '—'

# ── functional-marker completeness by category (this pipeline) ───────────────
func = load_tsv(os.path.join(DATA, 'novel_functional_ko_long.tsv'))
CATS = []
for r in func:
    if r['category'] not in CATS:
        CATS.append(r['category'])
cat_total = {c: len({r['gene'] for r in func if r['category'] == c}) for c in CATS}
cat_present = {m: {c: 0 for c in CATS} for m in NOVEL}
for r in func:
    if int(r['copies']) > 0:
        cat_present[r['MAG']][r['category']] += 1

def trend(vals):
    return '↑' if vals[-1] > vals[0] * 1.5 else ('↓' if vals[-1] < vals[0] / 1.5 else '→')

# ── assemble rows ────────────────────────────────────────────────────────────
rows = []
for m in NOVEL:
    s = S1[m]
    ab = [float(s[13 + i]) * 100 for i in range(5)]   # fraction -> %
    rows.append({
        'iterativeID': s[1], 'bin': m, 'phylum': s[2], 'class': s[3],
        'order': s[4], 'family': s[5], 'genus': s[6],
        'completeness_%': f'{s[7]:.2f}', 'contamination_%': f'{s[8]:.2f}',
        'size_Mb': f'{s[9] / 1e6:.2f}', 'contigs': s[10],
        'N50_kb': f'{s[12] / 1e3:.0f}', 'GC_%': f'{s[11] * 100:.1f}',
        'RED': f'{float(red[m]):.3f}' if red.get(m) else '',
        'novelty': 'species-novel (s__; no ANI hit; topology)',
        'method': method.get(m, ''),
        **{f'abund_{SAMPLES[i]}_%': f'{ab[i]:.3f}' for i in range(5)},
        'abund_trend': trend(ab),
        'N_role': role.get(m, ''),
        'denitrification': denit_str(m),
        'nosZ': nosz_str(m),
        **{f'func_{c}': f'{cat_present[m][c]}/{cat_total[c]}' for c in CATS},
    })

# ── write full TSV ───────────────────────────────────────────────────────────
os.makedirs(DATA, exist_ok=True)
cols = list(rows[0].keys())
with open(os.path.join(DATA, 'novel_bacteria_summary.tsv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=cols, delimiter='\t')
    w.writeheader()
    w.writerows(rows)

# ── write a readable markdown (two focused tables) ───────────────────────────
def md_table(headers, body):
    out = ['| ' + ' | '.join(headers) + ' |',
           '|' + '|'.join(['---'] * len(headers)) + '|']
    for r in body:
        out.append('| ' + ' | '.join(str(x) for x in r) + ' |')
    return '\n'.join(out)

with open(os.path.join(DATA, 'novel_bacteria_summary.md'), 'w') as f:
    f.write('# Species-novel MAGs of the CAN system — summary\n\n')
    f.write('Seven dereplicated MAGs were classified by GTDB-Tk to a named or '
            'placeholder genus but carry **no species assignment** (`s__`), with '
            'placement *fully defined by tree topology* and no ANI hit to any GTDB '
            'reference genome — i.e. each represents a candidate novel species.\n\n')

    f.write('### Table N1. Taxonomy, novelty and genome quality\n\n')
    f.write(md_table(
        ['iterativeID', 'bin', 'phylum', 'genus', 'compl %', 'contam %',
         'size (Mb)', 'contigs', 'N50 (kb)', 'GC %', 'RED'],
        [[r['iterativeID'], r['bin'], r['phylum'], r['genus'], r['completeness_%'],
          r['contamination_%'], r['size_Mb'], r['contigs'], r['N50_kb'],
          r['GC_%'], r['RED']] for r in rows]) + '\n\n')

    f.write('### Table N2. Abundance trajectory and N-cycle inventory\n\n')
    f.write(md_table(
        ['iterativeID', 'CAN_1 %', 'CAN_2 %', 'CAN_3 %', 'CAN_4 %', 'CAN_5 %',
         'trend', 'N-role', 'denitrification / DNRA', 'nosZ'],
        [[r['iterativeID'], r['abund_CAN_1_%'], r['abund_CAN_2_%'], r['abund_CAN_3_%'],
          r['abund_CAN_4_%'], r['abund_CAN_5_%'], r['abund_trend'], r['N_role'],
          r['denitrification'], r['nosZ']] for r in rows]) + '\n\n')

    f.write('### Table N3. Functional-marker completeness (genes present / genes screened)\n\n')
    f.write(md_table(
        ['iterativeID'] + CATS,
        [[r['iterativeID']] + [r[f'func_{c}'] for c in CATS] for r in rows]) + '\n\n')
    f.write('_Functional markers: KOfam (pyhmmer, adaptive thresholds) over Bakta '
            'proteomes. Genomic potential only — no metatranscriptome. CAZy/MEROPS/'
            'eggNOG deep annotation was not run, so secreted-protease and EPS '
            'inventories are out of scope; glycoside-hydrolase counts are KO-level only._\n')

# console preview
print('=== Table N1 preview ===')
for r in rows:
    print(f"  {r['iterativeID']:24s} {r['phylum']:16s} compl={r['completeness_%']}% "
          f"contam={r['contamination_%']}% {r['size_Mb']}Mb/{r['contigs']}ctg")
print('\n=== denitrification / abundance ===')
for r in rows:
    print(f"  {r['iterativeID']:24s} {r['abund_trend']} role={r['N_role']:16s} "
          f"denit=[{r['denitrification']}]  nosZ={r['nosZ']}")
print(f"\nwrote novel_bacteria_summary.tsv / .md to {DATA}")
