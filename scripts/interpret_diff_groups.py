from snakemake.script import snakemake
from archimedes.functions.ui_complements import subsetDF_index_targets
from archimedes.functions.dataflow import input_json, input_path, param, output_path
from archimedes.functions.manipulations import unique
import pandas as pd

### Tweakable Parameters
group_def_1 = input_json('group_def_1')
group_def_2 = input_json('group_def_1')
sample_metadata_file = input_path('pseudo_metadata')
sample_col = param('sample_id_column')
ct_col = param('cell_type_column')

cell_metadata = pd.read_csv(sample_metadata_file, index_col = 0)

### --- Quirk for this data --- ###
# Oops... re-making from the pseudobulk's names
import re
def empty_or_match(search):
    if search is None:
        return ''
    else:
        return search.group(1)
cell_metadata['manual.celltype'] = [empty_or_match(re.search(r"^.*__(.+)$", i)) for i in cell_metadata.index]
### --- End Quirk --- ###

# Gather samples in comparator groups
group_1_cells = subsetDF_index_targets(cell_metadata, group_def_1).index
group_2_cells = subsetDF_index_targets(cell_metadata, group_def_2).index

# Output indexes
group_1_cells.to_series().to_csv(output_path('group_1_inds'), header=False, index=False)
group_2_cells.to_series().to_csv(output_path('group_2_inds'), header=False, index=False)

# Output full metadata
group_1_meta = cell_metadata.loc[group_1_cells,:]
group_1_meta.to_csv(output_path('group_1_metadata'))
group_2_meta = cell_metadata.loc[group_2_cells,:]
group_2_meta.to_csv(output_path('group_2_metadata'))

# Summarize group contents
def to_md(text: str, indent: int = 0, extra_lines: int = 2):
    return '    '*indent + text + '\n'*extra_lines
def summarize_group(gnum: int, meta: pd.DataFrame):
    summary = to_md(f'### Group {gnum}:')
    npseudo = len(meta.index)
    samps = unique([str(x) for x in meta[sample_col]])
    cnums = list(meta['pseudobulk_cell_num'])
    cnums_low = [i for i in cnums if i < 100]
    cts = unique([str(x) for x in meta[ct_col]])
    cts_str = f': {", ".join(cts) if len(cts)<=5 else ""}'
    if npseudo < 1:
        summary += to_md(f'**ERRORS EXPECTED**: 0 pseudobulks met the defined criteria')
    if npseudo == 1:
        summary += to_md(f'**ERRORS EXPECTED**: Only 1 pseudobulk met the defined criteria')
    summary += to_md(f'This grouping contains {npseudo} pseduobulks from:')
    summary += to_md(f'- {len(samps)} sample(s)', 0, 1)
    summary += to_md(f'- {len(cts)} cell type(s){cts_str}')
    if npseudo <= 10 or len(samps) <= 10:
        summary += to_md(f'Note: Fewer than 10 {"pseudobulks, so statistical comparison has low N." if npseudo <= 10 else "original samples"}')
    if len(cnums_low)>0:
        summary += to_md(f'Note: {len(cnums_low)} pseudobulk(s) were built from fewer than 100 cells: {", ".join([str(x) for x in cnums_low])}')
    return summary
full_summary = summarize_group(1, group_1_meta) + summarize_group(2, group_2_meta)
with(open(snakemake.output['groups_summary'], 'w')) as file:
    file.write(full_summary)
