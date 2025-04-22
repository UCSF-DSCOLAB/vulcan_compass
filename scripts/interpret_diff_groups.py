from snakemake.script import snakemake
from archimedes.functions.ui_complements import subsetDF_index_targets
from archimedes.functions.dataflow import input_json, input_path, output_path
import pandas as pd

### Tweakable Parameters
group_def_1 = input_json('group_def_1')
group_def_2 = input_json('group_def_1')
sample_metadata_file = input_path('pseudo_metadata')

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

# Output
group_1_cells.to_series().to_csv(output_path('group_1_inds'), header=False, index=False)
group_2_cells.to_series().to_csv(output_path('group_2_inds'), header=False, index=False)
