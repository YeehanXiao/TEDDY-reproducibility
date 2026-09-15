#!/usr/bin/env python3
"""Summarize frozen FREDY-adapted calls, retaining novel transcript mappings."""
import argparse
import itertools
import re
import subprocess
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

REPS = ('rep1', 'rep2', 'rep3', 'rep4')
LOCI = ('Nelfa', 'Zfp352', 'Lmx1a', 'Cdk2ap1', 'Snai1', 'Pou6f2', 'Fam172a')
ATTR = re.compile(r'(\S+)\s+"([^"]*)"')


def records(path, minimum=1):
    with path.open() as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip('\r\n').split('\t')
            if len(fields) < minimum:
                raise ValueError(f'{path}:{number}: invalid column count')
            yield fields


def stable(value):
    return re.sub(r'\.\d+$', '', value) if value.startswith('ENS') else value


def strand(fields):
    values = set(fields[5:-1]) & {'+', '-'}
    if len(values) > 1:
        raise ValueError('Conflicting strand fields')
    return next(iter(values), '.')


def models(path):
    result = {}
    for f in records(path, 6):
        if f[2] not in ('transcript', 'mRNA', 'exon'):
            continue
        a = dict(ATTR.findall(f[-1]))
        tid = a.get('transcript_id')
        if not tid:
            raise ValueError(f'{path}: missing transcript_id')
        m = result.setdefault(tid, {'attrs': defaultdict(set), 'exons': set()})
        for key in ('gene_name', 'ref_gene_id', 'gene_id'):
            if a.get(key):
                m['attrs'][key].add(a[key])
        if f[2] == 'exon':
            start, end = int(f[3]), int(f[4])
            if start < 1 or end < start:
                raise ValueError(f'{path}: invalid exon coordinates for {tid}')
            m['exons'].add((f[0], start, end, strand(f)))
    if not result:
        raise ValueError(f'No transcript models: {path}')
    return result


def reference(path, bed):
    genes, transcripts = {}, {}
    by_name = defaultdict(set)
    with bed.open('w') as out:
        for f in records(path, 9):
            a = dict(ATTR.findall(f[8]))
            if 'gene_id' not in a:
                continue
            gid = stable(a['gene_id'])
            name = a.get('gene_name', gid)
            genes[gid] = name
            by_name[name].add(gid)
            if a.get('transcript_id'):
                tid = a['transcript_id']
                if stable(tid) in transcripts and transcripts[stable(tid)] != gid:
                    raise ValueError(f'Ambiguous reference transcript: {tid}')
                transcripts[stable(tid)] = gid
                if f[2] == 'exon' and int(f[4]) > int(f[3]):
                    # Native FREDY gene association uses GTF coordinates as BED fields.
                    out.write(f'{f[0]}\t{f[3]}\t{f[4]}\t{tid}\t0\t{f[6]}\n')
    if not genes or not transcripts or not bed.stat().st_size:
        raise ValueError('Reference annotation is empty or invalid')
    return genes, transcripts, by_name


def te_bed(path, bed):
    count = 0
    with bed.open('w') as out:
        for f in records(path, 9):
            start, end = int(f[3]) - 1, int(f[4])
            a, orientation = dict(ATTR.findall(f[8])), strand(f)
            name = next((a[k] for k in ('gene_name', 'repName', 'repeat_name',
                        'transcript_id', 'gene_id', 'family_id') if a.get(k)), None)
            name = name or f[8].strip().strip('"').strip(';').split(';')[0].strip().strip('"')
            if orientation not in ('+', '-') or not name or name == '.':
                continue
            if start < 0 or end <= start:
                raise ValueError(f'Invalid TE interval: {f[:5]}')
            out.write(f'{f[0]}\t{start}\t{end}\t{name}\t0\t{orientation}\n')
            count += 1
    if not count:
        raise ValueError('No stranded TE records were parsed')
    print(f'TE records: {count}', flush=True)


def intersect(a, b, bedtools, same_strand=False):
    command = [bedtools, 'intersect', '-a', str(a), '-b', str(b), '-wo']
    if same_strand:
        command.append('-s')
    with tempfile.TemporaryFile(mode='w+') as errors:
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=errors, text=True)
        try:
            for line in process.stdout:
                yield line.rstrip('\n').split('\t')
        finally:
            process.stdout.close()
            status = process.wait()
        if status:
            errors.seek(0)
            raise RuntimeError(f'bedtools failed ({status}): {errors.read()}')


def assign_genes(adapted, default, shared, genes, ref_tx, by_name, refbed, tmp, bedtools):
    mapping, source, pending = {}, {}, []
    for tid, m in adapted.items():
        a = m['attrs']
        candidates = {stable(g) for k in ('ref_gene_id', 'gene_id') for g in a[k]
                      if stable(g) in genes}
        if stable(tid) in ref_tx:
            candidates.add(ref_tx[stable(tid)])
        if not candidates:
            candidates = set().union(*(by_name.get(n, set()) for n in a['gene_name']))
        if tid in shared:
            if tid not in default or m['exons'] != default[tid]['exons']:
                raise ValueError(f'Cannot transfer default mapping: structure mismatch for {tid}')
            if candidates and len(shared[tid]) == 1 and candidates != shared[tid]:
                raise ValueError(f'Conflicting annotation/default mapping for {tid}')
            mapping[tid], source[tid] = shared[tid], 'default_most_shared'
        elif candidates:
            mapping[tid], source[tid] = candidates, 'reference_annotation'
        else:
            pending.append(tid)

    # Associate ALL remaining models, not only the manuscript loci.
    query = tmp / 'unmapped_exons.bed'
    with query.open('w') as out:
        for tid in pending:
            for chrom, start, end, orientation in sorted(adapted[tid]['exons']):
                if end > start:
                    out.write(f'{chrom}\t{start}\t{end}\t{tid}\t0\t{orientation}\n')
    scores = defaultdict(Counter)
    if query.stat().st_size:
        for f in intersect(query, refbed, bedtools):
            scores[f[3]][f[9]] += 1
    for tid in pending:
        counts = scores.get(tid, {})
        best = max(counts.values(), default=0)
        mapping[tid] = {ref_tx[stable(t)] for t, score in counts.items()
                        if score == best and best >= 2}
        source[tid] = 'max_shared_exon_pairs' if mapping[tid] else 'unresolved'
    return mapping, source


def first_exon(m):
    exons = m['exons']
    if not exons:
        return None
    contexts = {(e[0], e[3]) for e in exons}
    if len(contexts) != 1:
        raise ValueError('A transcript contains inconsistent chromosome/strand assignments')
    orientation = next(iter(contexts))[1]
    if orientation == '+':
        return min(exons, key=lambda e: (e[1], e[2]))
    if orientation == '-':
        return max(exons, key=lambda e: (e[2], e[1]))
    return None


def save(frame, out, name):
    frame.to_csv(out / name, sep='\t', index=False)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--default-root', type=Path, required=True)
    p.add_argument('--annotation', type=Path, required=True)
    p.add_argument('--te-gtf', type=Path, required=True)
    p.add_argument('--output-dir', type=Path, required=True)
    p.add_argument('--baseline-results', type=Path)
    p.add_argument('--bedtools', default='bedtools')
    args = p.parse_args()
    out = args.output_dir
    if args.baseline_results and out.resolve() == args.baseline_results.resolve():
        raise ValueError('Use a separate output directory for verification')
    out.mkdir(parents=True, exist_ok=True)
    hits, audit, qc = [], [], []
    checks = {('rep1', 'Lmx1a', 'MT2_Mm'), ('rep2', 'Snai1', 'MT2_Mm')}
    recovered_checks = set()
    with tempfile.TemporaryDirectory() as directory:
        tmp = Path(directory)
        genes, ref_tx, by_name = reference(args.annotation, tmp / 'reference.bed')
        te_bed(args.te_gtf, tmp / 'te.bed')
        for rep in REPS:
            sample = f'2cell{rep}'
            d = args.default_root / sample / 'chimeric'
            r = args.root / sample / 'chimeric'
            default, adapted = models(d / 'protein.gtf'), models(r / 'protein.gtf')
            initial = {f[0] for f in records(d / 'info.tsv', 5) if f[4] == 'Novel Initial'}
            missing = initial - adapted.keys()
            changed = {t for t in initial & adapted.keys()
                       if t not in default or default[t]['exons'] != adapted[t]['exons']}
            if missing or changed:
                raise ValueError(f'{rep}: default initial models missing={len(missing)}, changed={len(changed)}')
            shared = defaultdict(set)
            for f in records(d / 'most_shared.tsv', 3):
                gid = stable(f[1])
                if gid not in genes:
                    raise ValueError(f'Unknown reference gene: {gid}')
                shared[f[0]].add(gid)
            mapping, source = assign_genes(adapted, default, shared, genes, ref_tx,
                                         by_name, tmp / 'reference.bed', tmp, args.bedtools)
            selected, first = {}, {}
            with (tmp / 'first.bed').open('w') as handle:
                for tid, model in adapted.items():
                    candidates = mapping[tid]
                    if len(candidates) == 1:
                        gid = next(iter(candidates))
                        selected[tid] = (gid, genes[gid])
                    e = first_exon(model)
                    if e:
                        first[tid] = e
                        handle.write(f'{e[0]}\t{e[1]-1}\t{e[2]}\t{tid}\t0\t{e[3]}\n')
                    audit.append({'replicate': rep, 'transcript_id': tid,
                                  'gene_ids': ';'.join(sorted(candidates)),
                                  'gene_names': ';'.join(sorted({genes[g] for g in candidates})),
                                  'mapping_source': source[tid],
                                  'mapping_status': 'unique' if len(candidates) == 1 else
                                                    ('ambiguous' if candidates else 'unresolved')})
            hit_ids = set()
            for f in intersect(tmp / 'first.bed', tmp / 'te.bed', args.bedtools, True):
                overlap, tid = int(f[-1]), f[3]
                if overlap < 5:
                    continue
                hit_ids.add(tid)
                gid, name = selected.get(tid, ('', ''))
                hits.append({'replicate': rep, 'transcript_id': tid, 'gene_id': gid,
                             'gene_name': name, 'TE_name': f[9], 'chrom': f[0],
                             'first_exon_start': int(f[1])+1, 'first_exon_end': int(f[2]),
                             'strand': f[5], 'TE_start': int(f[7])+1, 'TE_end': int(f[8]),
                             'overlap_bp': overlap, 'mapping_source': source[tid]})
                if tid in initial and len(shared.get(tid, set())) == 1:
                    if mapping[tid] != shared[tid]:
                        raise ValueError(f'{rep}/{tid}: native unique gene association was lost')
                    recovered_checks.add((rep, name, f[9]))
            row = {'replicate': rep, 'protein_transcripts': len(adapted),
                   'default_initial_transcripts': len(initial), 'default_missing': len(missing),
                   'default_structure_changed': len(changed), 'gene_mapped_transcripts': len(selected),
                   'ambiguous_transcripts': sum(len(g) > 1 for g in mapping.values()),
                   'unresolved_transcripts': sum(not g for g in mapping.values()),
                   'transcripts_without_stranded_first_exon': len(adapted)-len(first),
                   'first_exon_TE_transcripts': len(hit_ids),
                   'gene_mapped_first_exon_TE_transcripts': len(hit_ids & selected.keys()),
                   'unresolved_or_ambiguous_first_exon_TE_transcripts': len(hit_ids-selected.keys())}
            qc.append(row)
            print(f'{rep}: {len(initial)} default initial models retained with identical exons; '
                  f'{len(selected)}/{len(adapted)} transcripts have unique gene mappings', flush=True)
    if checks - recovered_checks:
        raise ValueError(f'Regression failed: {sorted(checks-recovered_checks)}')
    events = pd.DataFrame(hits).drop_duplicates()
    called = events[events['gene_name'] != ''].copy()
    if called.empty:
        raise ValueError('No gene-mapped first-exon TE calls')
    sets = {rep: set(called.loc[called['replicate'] == rep, 'gene_name']) for rep in REPS}
    family = called.groupby('gene_name')['TE_name'].agg(lambda v: ';'.join(sorted(set(v))))
    recurrence = pd.DataFrame([{'gene_name': gene,
        **{f'{rep}_supported': gene in sets[rep] for rep in REPS},
        'TE_families': family[gene], 'n_replicates': sum(gene in sets[rep] for rep in REPS)}
        for gene in sorted(set().union(*sets.values()))])
    if args.baseline_results:
        baseline = pd.read_csv(args.baseline_results / 'fredy_adapted_gene_recurrence.tsv', sep='\t', dtype=str)
        for rep in REPS:
            old = set(baseline.loc[baseline[f'{rep}_supported'].str.lower().isin(['true', '1', 'yes']), 'gene_name'])
            if old - sets[rep]:
                raise ValueError(f'{rep}: previously counted gene calls lost: {sorted(old-sets[rep])[:20]}')
    summary = pd.DataFrame([{**{f'{r}_genes': len(sets[r]) for r in REPS},
        'union_genes': len(recurrence),
        'at_least_2_reps': int((recurrence.n_replicates >= 2).sum()),
        'at_least_3_reps': int((recurrence.n_replicates >= 3).sum()),
        'all_4_reps': int((recurrence.n_replicates == 4).sum())}])
    pairwise = pd.DataFrame([{'replicate_1': a, 'replicate_2': b,
        'intersection': len(sets[a] & sets[b]), 'union': len(sets[a] | sets[b]),
        'jaccard': len(sets[a] & sets[b])/len(sets[a] | sets[b])}
        for a, b in itertools.combinations(REPS, 2)])
    loci = []
    for gene in LOCI:
        row = {'gene_name': gene}
        for rep in REPS:
            families = set(called.loc[(called.gene_name == gene) & (called.replicate == rep), 'TE_name'])
            row[f'{rep}_supported'] = int(bool(families))
            row[f'{rep}_TE_families'] = ';'.join(sorted(families))
        row['n_replicates'] = sum(row[f'{r}_supported'] for r in REPS)
        row['TE_families'] = family.get(gene, '')
        loci.append(row)
    for frame, name in [(events, 'first_exon_events'), (pd.DataFrame(audit), 'transcript_gene_map'),
                        (pd.DataFrame(qc), 'mapping_qc'), (recurrence, 'gene_recurrence'),
                        (summary, 'overlap_summary'), (pairwise, 'pairwise_overlap'),
                        (pd.DataFrame(loci), 'locus_recovery')]:
        save(frame, out, f'fredy_adapted_{name}.tsv')
    print('\n', summary.to_string(index=False), '\n', pd.DataFrame(loci).to_string(index=False), sep='')
    print('\nFREDY_MAPPING_CHECKS_PASSED', flush=True)


if __name__ == '__main__':
    main()
