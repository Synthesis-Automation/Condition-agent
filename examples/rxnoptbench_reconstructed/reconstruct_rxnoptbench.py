#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Offline, source-preserving reconstruction of RxnOptBench MC_all JSONL.

Usage:
  python reconstruct_rxnoptbench.py --input multiple_choice_all_varying.jsonl --out results
Optional RDKit supplies syntax validation and canonical SMILES; nothing is fetched.
The raw source is preserved; source normalization is never treated as paper verification.
"""
from __future__ import annotations
import argparse
import collections
import copy
import csv
import hashlib
import json
import math
import re
import shutil
import sys
from pathlib import Path
from typing import Any
from urllib.parse import quote

try:
    from rdkit import Chem, rdBase
    rdBase.DisableLog('rdApp.error')
    rdBase.DisableLog('rdApp.warning')
    RDKIT_VERSION = rdBase.rdkitVersion
except ImportError:
    Chem = None
    RDKIT_VERSION = None

METRICS = ('yield', 'ee', 'dr', 'er', 'rr')
ALIGNED = ('yields', 'scores', 'effective_scores', 'option_relative_scores', 'option_source_entries')
SOURCE_URL = 'https://huggingface.co/datasets/songjhPKU/RxnOptBench/blob/main/data/jsonl/multiple_choice_all_varying.jsonl'
PAPER_URL = 'https://arxiv.org/html/2610.02242v1'
FRIENDLY = {
    'catalyst': 'catalysts', 'ligand': 'ligand', 'reagents': 'reagents',
    'solvents': 'solvents', 'temperature': 'reaction_temperature',
    'time': 'reaction_time', 'atmosphere': 'atmosphere',
    'light_condition': 'light_condition', 'electrolyte': 'Electrolyte',
    'anode': 'Anode(+)', 'cathode': 'Cathode(-)', 'current': 'Constant_Current',
    'current_density': 'Current_Density', 'cell_voltage': 'Constant_Cell_Voltage',
    'potential': 'Constant_Potential', 'charge_quantity': 'Constant_Quantity',
    'mode': 'Main_Mode', 'reaction_type_source': 'Reaction_Type',
    'stirring_speed': 'speed', 'pressure': 'pressure', 'pH': 'PH',
}


def dumps(x: Any) -> str:
    return json.dumps(x, ensure_ascii=False, sort_keys=True, separators=(',', ':'), allow_nan=False)


def digest(x: Any) -> str:
    return hashlib.sha256(dumps(x).encode('utf-8')).hexdigest()


def read_jsonl(path: Path) -> list[dict]:
    out = []
    with path.open(encoding='utf-8-sig') as f:
        for line_no, line in enumerate(f, 1):
            if not line.strip():
                continue
            try:
                obj = json.loads(line)
            except json.JSONDecodeError as exc:
                raise ValueError(f'{path.name}, line {line_no}: {exc}') from exc
            if not isinstance(obj, dict):
                raise TypeError(f'Line {line_no} is not an object')
            out.append(obj)
    if not out:
        raise ValueError('Input contains no records')
    return out


def write_json(path: Path, data: Any) -> None:
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def write_jsonl(path: Path, rows: list[dict]) -> None:
    with path.open('w', encoding='utf-8', newline='\n') as f:
        for row in rows:
            f.write(json.dumps(row, ensure_ascii=False, separators=(',', ':'), allow_nan=False) + '\n')


def write_csv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    fields = fields or list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w', encoding='utf-8-sig', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction='raise')
        writer.writeheader()
        for row in rows:
            writer.writerow({k: dumps(v) if isinstance(v, (dict, list)) else v for k, v in row.items()})


def normalized_key(key: str) -> str:
    return re.sub(r'[^a-z0-9]+', '_', key.lower().replace('+', ' plus ').replace('-', ' minus ')).strip('_')


def names(items: list[dict], key: str = 'content', sep: str = ' | ') -> str:
    # Human-readable convenience only. Full item alignment, empty strings and lists
    # remain in JSON and the long-format condition_items file.
    return sep.join(str(x.get(key, '')) for x in items if str(x.get(key, '')).strip())


def decode_doi(source_id: str) -> tuple[str | None, str]:
    if re.fullmatch(r'10\.\d{4,9}/\S+', source_id):
        return source_id.lower(), 'literal_doi_in_source'
    m = re.fullmatch(r'10_(\d{4,9})_(.+)', source_id)
    if not m:
        return None, 'unrecognized_source_identifier'
    registrant, suffix = m.groups()
    # Publisher-specific decoding is explicit inference, NOT resolver verification.
    if registrant == '1038' and re.fullmatch(r's\d+_\d{3}_\d+_[a-z0-9]+', suffix, re.I):
        return f'10.{registrant}/{suffix.replace("_", "-")}'.lower(), 'nature_hyphen_suffix_inferred'
    if registrant in {'1002', '1016', '1021', '1126'}:
        return f'10.{registrant}/{suffix.replace("_", ".")}'.lower(), 'period_separated_suffix_inferred'
    if registrant == '1039' and re.fullmatch(r'[a-z0-9]+', suffix, re.I):
        return f'10.{registrant}/{suffix}'.lower(), 'rsc_literal_suffix_inferred'
    return None, 'unrecognized_publisher_pattern'


def source_location(parent: str, source_id: str) -> dict:
    remainder = parent[len(source_id):] if parent.startswith(source_id) else ''
    m = re.fullmatch(r'_(?:(si|main)_)?(\d+)_table_(\d+)', remainder)
    return {
        'source_scope_token': (m.group(1) or 'unspecified') if m else 'unparsed',
        'source_page_index': int(m.group(2)) if m else None,
        'source_table_index': int(m.group(3)) if m else None,
    }


SMILES_CACHE: dict[str, dict] = {}

def smiles_info(s: str) -> dict:
    if s in SMILES_CACHE:
        return SMILES_CACHE[s]
    if not s.strip():
        v = {'status': 'missing', 'canonical_smiles': None, 'has_dummy_atom': None}
    elif Chem is None:
        v = {'status': 'not_checked', 'canonical_smiles': None, 'has_dummy_atom': None}
    else:
        try:
            mol = Chem.MolFromSmiles(s)
            if mol is None:
                v = {'status': 'invalid', 'canonical_smiles': None, 'has_dummy_atom': None}
            else:
                v = {'status': 'valid', 'canonical_smiles': Chem.MolToSmiles(mol, isomericSmiles=True),
                     'has_dummy_atom': any(a.GetAtomicNum() == 0 for a in mol.GetAtoms())}
        except Exception:
            v = {'status': 'invalid', 'canonical_smiles': None, 'has_dummy_atom': None}
    SMILES_CACHE[s] = v
    return v


def side_info(items: list[dict]) -> dict:
    info = [smiles_info(x['SMILES']) for x in items]
    complete = bool(items) and all(x['status'] != 'missing' for x in info)
    valid = complete and all(x['status'] == 'valid' for x in info)
    return {
        'raw': '.'.join(x['SMILES'] for x in items),
        'canonical': '.'.join(x['canonical_smiles'] for x in info) if valid else None,
        'complete': complete, 'valid': valid if Chem else None,
        'has_dummy': any(x['has_dummy_atom'] for x in info),
        'canonical_sorted': sorted(x['canonical_smiles'] for x in info) if valid else None,
    }


def numeric_temperature(text: str) -> tuple[float | None, str]:
    if not text:
        return None, 'missing_or_empty'
    t = text.strip().replace('−', '-').replace('–', '-').replace('℃', '°C').replace('º', '°')
    m = re.fullmatch(r'([+-]?\d+(?:\.\d+)?)\s*°\s*C', t, re.I)
    if m:
        return float(m.group(1)), 'explicit_single_celsius_value'
    if re.fullmatch(r'(?:r\.?\s*t\.?|room temperature|ambient(?: conditions)?)(?:\s*\((?:room temperature|r\.?\s*t\.?)\))?', t, re.I):
        return None, 'ambient_unspecified_numeric_value'
    return None, 'not_reduced_to_single_value'


def numeric_time(text: str) -> tuple[float | None, str]:
    if not text:
        return None, 'missing_or_empty'
    m = re.fullmatch(r'\s*(\d+(?:\.\d+)?)\s*(h|hr|hrs|hour|hours|min|mins|minute|minutes|s|sec|secs|second|seconds|d|day|days)\s*', text, re.I)
    if not m:
        return None, 'not_reduced_to_single_duration'
    number, unit = float(m.group(1)), m.group(2).lower()
    factor = 24 if unit in {'d', 'day', 'days'} else (1/60 if unit.startswith('min') else (1/3600 if unit.startswith('s') else 1))
    return number * factor, 'explicit_single_duration_converted_to_hours'


def reconstruct(source: Path, out_dir: Path) -> dict:
    records = read_jsonl(source)
    raw_sha = hashlib.sha256(source.read_bytes()).hexdigest()
    out_dir.mkdir(parents=True, exist_ok=True)
    raw_dir = out_dir / 'raw'
    raw_dir.mkdir(exist_ok=True)
    copy_path = raw_dir / source.name
    if source.resolve() != copy_path.resolve():
        shutil.copyfile(source, copy_path)
    issues: list[dict] = []
    rows: list[dict] = []
    tables: list[dict] = []
    condition_items: list[dict] = []
    molecule_occurrences: list[dict] = []
    structure_uses: dict[str, dict] = {}
    table_case_counts = collections.Counter(r['meta']['parent_table_id'].casefold() for r in records)
    source_ids_by_ref: dict[str, set] = collections.defaultdict(set)
    for rec in records:
        source_ids_by_ref[rec['meta']['doi'].casefold()].add(rec['meta']['doi'])
    condition_keys = sorted({k for r in records for c in [r['input']['constant_conditions'], *r['options']] for k in c})
    key_map = {k: normalized_key(k) for k in condition_keys}
    if len(set(key_map.values())) != len(key_map):
        raise ValueError('Condition column normalization collision')
    if len({r['id'] for r in records}) != len(records):
        raise ValueError('Duplicate benchmark question IDs; ambiguous provenance')
    if len({r['meta']['parent_table_id'] for r in records}) != len(records):
        raise ValueError('Duplicate exact table IDs; inspect input before reconstruction')

    def issue(level: str, scope: str, ident: str, code: str, detail: Any) -> None:
        issues.append({'severity': level, 'scope': scope, 'record_id': ident, 'issue_code': code, 'detail': detail})

    def track(s: str, name: str, role: str, table: str, row: str | None) -> None:
        if s not in structure_uses:
            structure_uses[s] = {'names': set(), 'roles': set(), 'tables': set(), 'rows': set(), 'occurrences': 0}
        u = structure_uses[s]
        if name: u['names'].add(name)
        u['roles'].add(role)
        u['tables'].add(table)
        if row: u['rows'].add(row)
        u['occurrences'] += 1

    for table_no, rec in enumerate(records, 1):
        table_ref = f'T{table_no:04d}'
        meta, inp, opts = rec['meta'], rec['input'], rec['options']
        parent, qid, sid = meta['parent_table_id'], rec['id'], meta['doi']
        n = len(opts)
        if rec['question_type'] != 'all_varying':
            raise ValueError(f'Unexpected question_type at {qid}')
        for key in ALIGNED:
            if len(meta[key]) != n:
                raise ValueError(f'{qid}: {key} length does not equal options length')
        if any(not isinstance(i, int) or i < 0 or i >= n for i in rec['answer']):
            raise ValueError(f'Invalid answer option at {qid}')
        if any(meta['scores'][i]['yield'] != meta['yields'][i] for i in range(n)):
            raise ValueError(f'Mismatched source yield arrays at {qid}')
        const = inp['constant_conditions']
        source_total = meta['total_entries_in_source_table']
        pointer_indices = [x['row_index'] for x in meta['option_source_entries']]
        if len(set(pointer_indices)) != len(pointer_indices) or any(not isinstance(i, int) or i < 0 or i >= source_total for i in pointer_indices):
            raise ValueError(f'Ambiguous or out-of-range source row pointers at {qid}')
        missing_indices = sorted(set(range(source_total)) - set(pointer_indices))
        doi, decode_rule = decode_doi(sid)
        refid = sid.casefold()
        location = source_location(parent, sid)
        a, b = side_info(inp['reactants']), side_info(inp['products'])
        reaction_raw = f"{a['raw']}>>{b['raw']}"
        reaction_can = f"{a['canonical']}>>{b['canonical']}" if a['canonical'] and b['canonical'] else None
        group_id = 'RX_' + digest([a['canonical_sorted'], b['canonical_sorted']])[:20] if reaction_can else None
        table_flags = []
        if missing_indices:
            issue('info', 'table', table_ref, 'source_entries_not_retained', {'count': len(missing_indices), 'row_indices': missing_indices})
        if table_case_counts[parent.casefold()] > 1:
            table_flags.append('case_colliding_table_identifier_retained')
            issue('warning', 'table', table_ref, 'case_colliding_table_identifier_retained', {'parent_table_id': parent, 'casefold_id': parent.casefold()})
        if not doi:
            issue('warning', 'table', table_ref, 'doi_not_decoded', sid)
        if not a['complete'] or not b['complete'] or a['valid'] is False or b['valid'] is False:
            table_flags.append('reaction_structure_needs_review')
            issue('warning', 'table', table_ref, 'reaction_structure_needs_review', reaction_raw)
        for side in ('reactants', 'products'):
            for j, item in enumerate(inp[side]):
                si = smiles_info(item['SMILES'])
                molecule_occurrences.append({'table_record_id': table_ref, 'parent_table_id': parent,
                    'reference_id': refid, 'side': side, 'molecule_index': j, 'name': item['content'],
                    'smiles_source': item['SMILES'], 'smiles_canonical': si['canonical_smiles'],
                    'smiles_status': si['status'], 'has_dummy_atom': si['has_dummy_atom']})
                track(item['SMILES'], item['content'], side, table_ref, None)
        common = {
            'table_record_id': table_ref, 'parent_table_id': parent, 'benchmark_question_id': qid,
            'source_record_number': table_no, 'reference_id': refid, 'doi_source_id': sid,
            'doi_inferred': doi, 'doi_url_inferred': 'https://doi.org/' + quote(doi, safe='/.-') if doi else None,
            'doi_decode_rule': decode_rule, 'doi_resolver_verified': False,
            'source_file': meta['source_file'], **location,
            'reaction_smiles': reaction_raw, 'reaction_smiles_canonical': reaction_can,
            'reactant_smiles': a['raw'], 'product_smiles': b['raw'],
            'reactant_smiles_canonical': a['canonical'], 'product_smiles_canonical': b['canonical'],
            'reaction_group_id': group_id,
            'reactant_names': names(inp['reactants']), 'product_names': names(inp['products']),
            'reaction_smiles_valid': (a['valid'] and b['valid']) if Chem else None,
            'reaction_has_dummy_atoms': a['has_dummy'] or b['has_dummy'],
            'main_product_smiles_source': meta['main_product_smiles'],
        }
        table_rows = []
        max_yield = max(meta['yields'])
        max_effective = max(meta['effective_scores'])
        for i, opt in enumerate(opts):
            if set(opt) != set(rec['varying_keys']):
                raise ValueError(f'Varying keys differ from option keys at {qid}, option {i}')
            if set(const).intersection(opt):
                raise ValueError(f'Constant/varying condition overlap at {qid}, option {i}')
            merged = {**copy.deepcopy(const), **copy.deepcopy(opt)}
            row_id = f'{table_ref}_O{i:04d}'
            pointer, score = meta['option_source_entries'][i], meta['scores'][i]
            if any(not math.isfinite(v) for v in score.values() if isinstance(v, (float, int))):
                raise ValueError(f'Non-finite outcome at {row_id}')
            flags = table_flags.copy()
            invalid_condition_items, missing_condition_smiles = [], 0
            for k, items in merged.items():
                for j, item in enumerate(items):
                    smi = item['SMILES']
                    si = smiles_info(smi)
                    track(smi, item['content'], f'condition:{k}', table_ref, row_id)
                    condition_items.append({
                        'screening_row_id': row_id, 'table_record_id': table_ref,
                        'condition_key': k, 'condition_origin': 'varying' if k in opt else 'constant',
                        'item_index': j, 'content_source': item['content'], 'smiles_source': smi,
                        'smiles_canonical': si['canonical_smiles'], 'smiles_status': si['status'],
                        'has_dummy_atom': si['has_dummy_atom'],
                    })
                    if si['status'] == 'invalid':
                        invalid_condition_items.append({'condition_key': k, 'item_index': j, 'content': item['content'], 'SMILES': smi})
                    if si['status'] == 'missing':
                        missing_condition_smiles += 1
                    if si['has_dummy_atom']:
                        flags.append('condition_smiles_contains_wildcard')
            if invalid_condition_items:
                flags.append('condition_smiles_parse_failure')
                issue('warning', 'row', row_id, 'condition_smiles_parse_failure', invalid_condition_items)
            if 'condition_smiles_contains_wildcard' in flags:
                issue('warning', 'row', row_id, 'condition_smiles_contains_wildcard', 'A condition structure contains a dummy/wildcard atom; source value retained.')
            row = {'screening_row_id': row_id, **common,
                'option_index': i, 'row_index': pointer['row_index'], 'entry_idx': pointer['entry_idx'], 'rxn_id': pointer['rxn_id']}
            for friendly, key in FRIENDLY.items():
                row[friendly] = names(merged.get(key, []))
                if friendly in {'catalyst', 'ligand', 'reagents', 'solvents'}:
                    row[friendly + '_smiles'] = names(merged.get(key, []), 'SMILES')
            row['temperature_C'], row['temperature_parse_status'] = numeric_temperature(row['temperature'])
            row['time_h'], row['time_parse_status'] = numeric_time(row['time'])
            row.update({k: score.get(k) for k in METRICS})
            row.update({
                'effective_score': meta['effective_scores'][i], 'relative_score': meta['option_relative_scores'][i],
                'is_benchmark_best': i in rec['answer'], 'is_max_yield_retained': score['yield'] == max_yield,
                'is_max_effective_score_retained': meta['effective_scores'][i] == max_effective,
                'is_global_best_entry_source': str(pointer['entry_idx']) in {str(x) for x in meta['global_best_entry_indices']},
                'varying_keys': rec['varying_keys'], 'condition_keys': sorted(merged),
                'empty_condition_keys': sorted(k for k,v in merged.items() if v == []),
                'missing_condition_keys': sorted(set(condition_keys) - set(merged)),
                'invalid_condition_smiles_count': len(invalid_condition_items),
                'empty_condition_smiles_item_count': missing_condition_smiles,
                'quality_flags': sorted(set(flags)),
                'conditions_json': merged, 'constant_conditions_json': const, 'varying_conditions_json': opt,
                'reactants_json': inp['reactants'], 'products_json': inp['products'],
                'scores_json': score, 'source_entry_json': pointer,
            })
            # All source condition fields get convenient flattened columns too.
            for k in condition_keys:
                row['cond__' + key_map[k]] = names(merged.get(k, []))
                row['cond_smiles__' + key_map[k]] = names(merged.get(k, []), 'SMILES')
            row['content_signature'] = digest({'reference_id': refid,
                'reactants': inp['reactants'], 'products': inp['products'], 'conditions': merged, 'scores': score})
            table_rows.append(row)
        # Present experimental source order, NOT the shuffled benchmark option order.
        table_rows.sort(key=lambda x: (x['row_index'], x['option_index']))
        rows.extend(table_rows)
        table = {**common,
            'n_reconstructed_options': n, 'total_entries_in_source_table': source_total,
            'n_source_entries_not_retained': source_total - n, 'missing_source_row_indices': missing_indices,
            'varying_keys': rec['varying_keys'], 'constant_condition_keys': sorted(const),
            'all_condition_keys': sorted(set(const) | set(rec['varying_keys'])),
            'yield_min_retained': min(meta['yields']), 'yield_max_retained': max_yield,
            'yield_mean_retained': sum(meta['yields']) / n,
            'best_yield_source': meta['best_yield'], 'best_effective_score_source': meta['best_effective_score'],
            'best_option_indices_source': rec['answer'], 'best_entry_indices_source': meta['best_entry_indices'],
            'global_best_yield_source': meta['global_best_yield'],
            'global_best_effective_score_source': meta['global_best_effective_score'],
            'global_best_entry_indices_source': meta['global_best_entry_indices'],
            'n_rows_with_invalid_condition_smiles': sum(x['invalid_condition_smiles_count'] > 0 for x in table_rows),
            'quality_flags': table_flags, 'constant_conditions_json': const,
            'reactants_json': inp['reactants'], 'products_json': inp['products'],
            'source_record_json': rec,
        }
        tables.append(table)

    signatures = collections.defaultdict(list)
    for r in rows: signatures[r['content_signature']].append(r['screening_row_id'])
    for r in rows:
        matches = signatures[r['content_signature']]
        r['identical_content_group_size'] = len(matches)
        if len(matches) > 1:
            r['quality_flags'] = sorted(set(r['quality_flags'] + ['identical_released_content_elsewhere']))
    for sig, matches in signatures.items():
        if len(matches) > 1:
            issue('info', 'group', sig, 'identical_released_content_retained', {'row_ids': matches, 'count': len(matches)})

    references = []
    for refid, source_ids in sorted(source_ids_by_ref.items()):
        ts = [t for t in tables if t['reference_id'] == refid]
        rs = [r for r in rows if r['reference_id'] == refid]
        first = ts[0]
        references.append({
            'reference_id': refid, 'doi_inferred': first['doi_inferred'],
            'doi_url_inferred': first['doi_url_inferred'], 'doi_decode_rule': first['doi_decode_rule'],
            'doi_resolver_verified': False, 'source_id_variants': sorted(source_ids),
            'n_source_id_variants': len(source_ids), 'n_table_records': len(ts),
            'n_reconstructed_rows': len(rs),
            'total_entries_in_source_table_sum': sum(t['total_entries_in_source_table'] for t in ts),
            'n_source_entries_not_retained': sum(t['n_source_entries_not_retained'] for t in ts),
            'table_record_ids': [t['table_record_id'] for t in ts],
            'parent_table_ids': [t['parent_table_id'] for t in ts],
            'bibliographic_metadata_status': 'Title, authors, journal and publication date are not supplied in the uploaded file.',
        })
    smiles_audit = []
    for s, use in sorted(structure_uses.items()):
        if not s.strip():
            continue
        si = smiles_info(s)
        smiles_audit.append({'smiles_source': s, 'smiles_canonical': si['canonical_smiles'],
            'smiles_status': si['status'], 'has_dummy_atom': si['has_dummy_atom'],
            'source_names': sorted(use['names']), 'source_roles': sorted(use['roles']),
            'n_table_records': len(use['tables']), 'n_screening_rows_as_condition': len(use['rows']),
            'n_occurrences_in_exported_molecule_and_condition_item_tables': use['occurrences'],
            'table_record_ids': sorted(use['tables']), 'screening_row_ids_as_condition': sorted(use['rows'])})

    # Exhaustive round-trip comparison to the supplied benchmark (NOT original PDFs).
    row_by_source = {(r['benchmark_question_id'], r['option_index']): r for r in rows}
    for rec in records:
        for i, opt in enumerate(rec['options']):
            r = row_by_source[(rec['id'], i)]
            assert r['conditions_json'] == {**rec['input']['constant_conditions'], **opt}
            assert r['varying_conditions_json'] == opt
            assert r['constant_conditions_json'] == rec['input']['constant_conditions']
            assert r['reactants_json'] == rec['input']['reactants']
            assert r['products_json'] == rec['input']['products']
            assert r['scores_json'] == rec['meta']['scores'][i]
            assert r['source_entry_json'] == rec['meta']['option_source_entries'][i]
            assert r['effective_score'] == rec['meta']['effective_scores'][i]
            assert r['relative_score'] == rec['meta']['option_relative_scores'][i]
            assert r['is_benchmark_best'] == (i in rec['answer'])
    assert [t['source_record_json'] for t in tables] == records
    assert len({r['screening_row_id'] for r in rows}) == len(rows)
    assert sum(t['n_reconstructed_options'] for t in tables) == len(rows)
    assert sum(t['n_source_entries_not_retained'] for t in tables) == sum(t['total_entries_in_source_table'] for t in tables) - len(rows)

    report = {
        'input_file': source.name, 'input_sha256': raw_sha, 'input_bytes': source.stat().st_size,
        'input_source_url_context_only': SOURCE_URL, 'paper_url_context_only': PAPER_URL,
        'external_sources_fetched': False, 'rdkit_version': RDKIT_VERSION,
        'n_table_records': len(tables), 'n_screening_rows': len(rows),
        'n_source_identifiers_exact': len({t['doi_source_id'] for t in tables}),
        'n_source_identifiers_case_insensitive': len(references),
        'n_table_identifiers_exact': len({t['parent_table_id'] for t in tables}),
        'n_table_identifiers_case_insensitive': len({t['parent_table_id'].casefold() for t in tables}),
        'source_identifier_case_variant_groups': [r['source_id_variants'] for r in references if r['n_source_id_variants'] > 1],
        'case_colliding_table_records': [t['table_record_id'] for t in tables if 'case_colliding_table_identifier_retained' in t['quality_flags']],
        'sum_source_table_entry_counts': sum(t['total_entries_in_source_table'] for t in tables),
        'n_source_entries_not_retained': sum(t['n_source_entries_not_retained'] for t in tables),
        'n_tables_with_unretained_source_entries': sum(t['n_source_entries_not_retained'] > 0 for t in tables),
        'fraction_of_source_entry_counts_retained': len(rows) / sum(t['total_entries_in_source_table'] for t in tables),
        'n_source_condition_keys': len(condition_keys), 'source_condition_keys': condition_keys,
        'n_condition_item_records': len(condition_items), 'n_reaction_molecule_occurrences': len(molecule_occurrences),
        'n_unique_source_reaction_molecule_smiles': len({m['smiles_source'] for m in molecule_occurrences}),
        'n_tables_with_valid_reaction_smiles': sum(t['reaction_smiles_valid'] is True for t in tables),
        'n_rows_with_valid_reaction_smiles': sum(r['reaction_smiles_valid'] is True for r in rows),
        'n_canonical_reaction_groups': len({t['reaction_group_id'] for t in tables if t['reaction_group_id']}),
        'n_unique_nonempty_condition_smiles': len({c['smiles_source'] for c in condition_items if c['smiles_source'].strip()}),
        'n_unique_invalid_nonempty_condition_smiles': len({c['smiles_source'] for c in condition_items if c['smiles_status'] == 'invalid'}),
        'n_invalid_condition_smiles_occurrences': sum(c['smiles_status'] == 'invalid' for c in condition_items),
        'n_rows_with_invalid_condition_smiles': sum(r['invalid_condition_smiles_count'] > 0 for r in rows),
        'n_tables_with_invalid_condition_smiles': sum(t['n_rows_with_invalid_condition_smiles'] > 0 for t in tables),
        'n_condition_items_with_wildcard': sum(c['has_dummy_atom'] is True for c in condition_items),
        'outcome_non_null_row_counts': {k: sum(r[k] is not None for r in rows) for k in METRICS},
        'n_zero_yield_rows': sum(r['yield'] == 0 for r in rows),
        'n_rows_yield_le_10': sum(r['yield'] <= 10 for r in rows),
        'n_rows_with_explicit_numeric_temperature_C': sum(r['temperature_C'] is not None for r in rows),
        'n_rows_with_explicit_numeric_time_h': sum(r['time_h'] is not None for r in rows),
        'n_identical_released_content_groups': sum(len(v) > 1 for v in signatures.values()),
        'n_extra_rows_with_identical_released_content': sum(len(v)-1 for v in signatures.values()),
        'n_rows_in_identical_released_content_groups': sum(len(v) for v in signatures.values() if len(v) > 1),
        'n_dois_decoded_but_not_resolver_verified': sum(bool(r['doi_inferred']) for r in references),
        'condition_key_coverage_rows': {k: {'key_present': sum(k in r['conditions_json'] for r in rows),
            'nonempty_list': sum(bool(r['conditions_json'].get(k)) for r in rows),
            'explicit_empty_list': sum(r['conditions_json'].get(k) == [] for r in rows)} for k in condition_keys},
        'quality_issue_counts': dict(collections.Counter(i['issue_code'] for i in issues)),
        'round_trip_validation': {'options_checked': len(rows), 'table_source_records_checked': len(tables),
            'all_reactants_products_conditions_scores_provenance_and_answers_match_input': True},
        'validation_not_performed': ['Original article / supporting information verification',
            'Public reference-subset cross-validation (not uploaded)', 'DOI resolver / Crossref verification',
            'Chemical identity, atom balance, atom mapping, or mechanistic validation'],
    }
    write_jsonl(out_dir / 'rxnoptbench_screening_rows.jsonl', rows)
    write_csv(out_dir / 'rxnoptbench_screening_rows.csv', rows)
    write_jsonl(out_dir / 'rxnoptbench_screening_tables.jsonl', tables)
    # Entire source record remains in table JSONL/raw, not the table-summary CSV.
    write_csv(out_dir / 'rxnoptbench_screening_tables.csv', [{k:v for k,v in t.items() if k != 'source_record_json'} for t in tables])
    write_csv(out_dir / 'rxnoptbench_references.csv', references)
    write_jsonl(out_dir / 'rxnoptbench_references.jsonl', references)
    write_csv(out_dir / 'rxnoptbench_condition_items.csv', condition_items)
    write_csv(out_dir / 'rxnoptbench_reaction_molecules.csv', molecule_occurrences)
    write_csv(out_dir / 'rxnoptbench_smiles_validation.csv', smiles_audit)
    write_csv(out_dir / 'rxnoptbench_quality_issues.csv', issues)
    write_json(out_dir / 'rxnoptbench_condition_key_map.json', key_map)
    write_json(out_dir / 'rxnoptbench_quality_report.json', report)
    # Round-trip exported representations as a final integrity check.
    assert read_jsonl(out_dir / 'rxnoptbench_screening_rows.jsonl') == rows
    with (out_dir / 'rxnoptbench_screening_rows.csv').open(encoding='utf-8-sig', newline='') as f:
        csv_rows = list(csv.DictReader(f))
    assert len(csv_rows) == len(rows)
    for parsed, original in zip(csv_rows, rows):
        assert json.loads(parsed['conditions_json']) == original['conditions_json']
        assert json.loads(parsed['scores_json']) == original['scores_json']
        assert parsed['reaction_smiles'] == original['reaction_smiles']
        assert parsed['entry_idx'] == original['entry_idx']
    return report


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', type=Path, required=True)
    p.add_argument('--out', type=Path, default=Path('rxnoptbench_reconstructed'))
    args = p.parse_args()
    if not args.input.is_file():
        p.error(f'Input file not found: {args.input}')
    result = reconstruct(args.input.resolve(), args.out.resolve())
    print(json.dumps(result, ensure_ascii=False, indent=2))

if __name__ == '__main__':
    main()
