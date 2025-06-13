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
    is_empty,
)


def _zero_pad(x, n):
    if len(x) > n:
        raise ValueError('list too long')
    return x + ([0]*(n-len(x)))


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
    simple_map('quantity_type', 'MT', 'required'),
    identity_map('comments', 'optional'),
    oneway_func_map('num_reaction_ids', lambda d: len(d['NT']), 'required', 'forward'),
    oneway_func_map(
        'reaction_ids', lambda d: _zero_pad(d['NT'], 5), 'required', 'forward'
    ),
    oneway_func_map(
        'NT',
        lambda d: d['reaction_ids'][:d['num_reaction_ids']],
        'required', 'backward'
    ),
    identity_map('NNCOX', 'required' ),
    simple_map('ENF', 'ENFF', 'optional'),
    identity_map('NENF', 'optional'),
    simple_map('EPA', 'EPAF', 'required'),
    identity_map('NETG', 'required'),
    simple_map('energies', 'E', 'required'),
    simple_map('measured_values', 'CSS', 'required'),
    simple_map('uncertainties', 'CO', 'required', transpose=True),
    simple_map('NCST', 'NCSST', 'optional'),
    identity_map('NEC', 'optional'),
    identity_map('FCFC', 'optional'),
]


def map_priorblock(priorblock: dict, direction: str='forward', do_reduce: bool=False):
    success = False
    for midx, mappings in enumerate((GMAPY_DATPY_PRIOR_MAPPINGS, GMAPY_DATPY_FISSION_MAPPINGS)):
        try:
            new_priorblock = map_dict(priorblock, mappings, direction)
            if midx == 0 and do_reduce and direction == 'backward':
                new_priorblock['EN'] = new_priorblock['EN'][1:-1]
                new_priorblock['CS'] = new_priorblock['CS'][1:-1]
            success = True
            break
        except Exception as Exc:
            pass
    if not success:
        raise ValueError('unable to map priorblock')
    return new_priorblock


def map_priorblocks(priorblocks: list, direction: str='forward', do_reduce: bool=False):
    new_priorblocks = []
    for idx, block in priorblocks.items():
        new_priorblock = map_priorblock(block, direction, do_reduce)
        new_priorblocks.append(new_priorblock)
    return new_priorblocks


def map_dataset(dataset: dict, direction: str='forward'):
    return map_dict(dataset, GMAPY_DATPY_DATASET_MAPPINGS, direction)


def map_datablocks(datablocks: list, direction: str='forward'):
    if direction not in ('forward', 'backward'):
        raise ValueError('direction must be `forward` or `backward`')

    new_datablocks = []
    for block in tuple(datablocks):

        datasets = block['datasets'] if direction == 'forward' else block

        if direction == 'forward':
            datasets = block['datasets']
            new_datasets = [map_dataset(ds, direction) for ds in datasets]
            # special-casing for moving the correlation matrix
            # from block-level into the last dataset of the block
            if not is_empty(block.get('ECOR')):
                new_datasets[-1]['cormat'] = block['ECOR']
            new_block = new_datasets

        elif direction == 'backward':
            datasets = block
            new_datasets = [map_dataset(ds, direction) for ds in datasets]
            new_block = {
                'type': 'legacy-experiment-datablock',
                'datasets': new_datasets
            }
            if not is_empty(datasets[-1].get('cormat')):
                new_block['ECOR'] = datasets[-1]['cormat']

        new_datablocks.append(new_block)

    return new_datablocks
