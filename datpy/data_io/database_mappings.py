from copy import deepcopy
import json
from ..datamodels.models import (
    Dataset,
    ReactionPrior,
    FissionSpectrum,
)
from .mapping_utils import (
    simple_augmentation,
    simple_map,
    identity_map,
    oneway_func_map,
    map_dict,
)
from .mapping_utils import map_dict


GMAPY_DATPY_FISSION_MAPPINGS = [
    simple_augmentation('type', 'legacy-fission-spectrum', 'backward'),
    simple_augmentation('label', 'spectrum', 'forward'),
    simple_map('energies', 'ENFIS', 'required'),
    simple_map('spectrum_values', 'FIS', 'required'),
]


GMAPY_DATPY_PRIOR_MAPPINGS = [
    simple_augmentation('type', 'legacy-prior-cross-section', 'backward'),
    simple_map('label', 'CLAB', 'required'),
    simple_map('reaction_id', 'ID', 'required'),
    simple_map('energies', 'EN', 'required'), 
    simple_map('cross_sections', 'CS', 'required'),
]


GMAPY_DATPY_DATASET_MAPPINGS = [
    simple_augmentation('type', 'legacy-experiment-dataset', 'backward'),
    simple_map('dataset_id', 'NS', 'required'),
    simple_map('year', 'YEAR', 'required'),
    simple_map('author', 'CLABL', 'required'),
    simple_map('pubref', 'BREF', 'required'),
    simple_map('tag', 'TAG', 'required'),
    simple_map('quantity_type', 'MT', 'required')
    identity_map('comments', 'optional'),
    oneway_func_map('num_reaction_ids', lambda d: len(d['NT']), 'required'),
    simple_map('reaction_ids', 'NT', 'required'),
    simple_map('ENF', 'ENFF', 'optional'),
    identity_map('NENF', 'optional'),
    simple_map('EPA', 'EPAF', 'required'),
    identity_map('NETG', 'required'),
    simple_map('energies', 'E', 'required'),
    simple_map('measured_values', 'CSS', 'required'),
    simple_map('uncertainties', 'CO', 'required'),
    simple_map('NCST', 'NCSST', 'required'),
    identity_map('NEC', 'required'),
    identity_map('FCFC', 'required'),
    simple_map('cormat', 'ECOR', 'required'),
]


def map_priorblock(priorblock: dict, direction: str='forward'):
    success = False
    for mappings in (GMAPY_DATPY_PRIOR_MAPPINGS, GMAPY_DATPY_FISSION_MAPPINGS):
        try:
            new_priorblock = map_dict(priorblock, mappings, direction)
            success = True
            break
        except:
            pass
    if not success:
        raise ValueError('unable to map priorblock')
    return new_priorblock


def map_priorblocks(priorblocks: list, direction: str='forward'):
    new_priorblocks = []
    for block in priorblocks:
        new_priorblock = map_priorblock(block, direction)
        new_priorblocks.append(new_priorblock)
    return new_priorblocks


def map_dataset(dataset: dict, direction: str='forward'):
    return map_dict(dataset, GMAPY_DATPY_DATASET_MAPPINGS, direction)


def map_datablocks(datablocks: list, direction: str='forward'):
    new_datablocks = deepcopy(datablocks) 
    for block in tuple(new_datablocks):
        datasets = block['datasets']
        # special-casing for moving the correlation matrix
        # from block-level into the last dataset of the block
        if direction == 'forward' and 'ECOR' in block:
            datasets[-1]['ECOR'] = block.pop('ECOR')
        # here the generic mapping
        new_datasets = [map_dataset(ds, direction) for ds in datasets]
        block['datasets'] = new_datasets
        # special-casing for moving the correlation matrix
        # from last dataset to block level
        if direction == 'backward' and 'ECOR' in new_datasets[-1]:
            block['ECOR'] = new_datasets[-1].pop('ECOR')
    return new_datablocks
