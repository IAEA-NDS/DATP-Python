import os
import io
import sys
import json
from pathlib import Path
from typing import Optional, Dict
from .data_io.legacy.input_output import (
    copy_gma_controls,
)
from .data_io.legacy.data_input import (
    read_apriori,
    read_datablocks,
)
from .data_io.legacy.input_output import (
    copy_gma_controls,
)
from .data_io.legacy..gma_output import (
    write_gmadb,
    write_prior,
)
from .reduction import reduce_datablocks
import argparse
import logging
from .datamodels.models import (
    Dataset,
    ReactionPrior,
)
from .reduction import reduce_datablocks


logger = logging.getLogger(__name__)


# extract input data and convert to modern dict/json

def extract_datasets(gmdata_crd_content: str) -> Dict[int, dict]:

    datasets = {}
    fileobj = io.StringIO(gmdata_crd_content)
    datablocks = read_datablocks(fileobj)
    for block in datablocks:
        for dataset in block:
            curid = dataset.dataset_id
            if curid in datasets:
                if len(block) == 1:
                    newid = curid + 1
                    while newid in datasets:
                        newid += 1
                    logger.warning(
                        f'Dataset {curid} already exists, use new id {newid}',
                    )
                    curid = newid
                    dataset.dataset_id = curid
                else:
                    raise KeyError(
                        f'Dataset {curid} already exists in block with '
                        'with several other datasets, reassignment not possible.'
                    )
            datasets[curid] = dataset.dict(use_arrays=False)
    return datasets


def extract_prior(prior_content: str) -> dict:
    fileobj = io.StringIO(prior_content)
    copy_gma_controls(fileobj, file_IO2=None, gma_file_handle=None)
    reaction_prior = read_apriori(fileobj)
    return reaction_prior.dict()


def extract_spectrum(prior_content: str) -> dict:
    fileobj = io.StringIO(prior_content)
    ret = copy_gma_controls(fileobj)
    return ret['spectrum']


# reduction functionality

def reduce_datasets(datasets: Dict[int, dict], prior: Dict[int, dict]) -> dict:
    reaction_prior = ReactionPrior(**prior)
    reduced_datasets = {}
    for dataset_id, raw_dataset in datasets.items():
        dataset = Dataset(**raw_dataset)
        if int(dataset_id) != dataset.dataset_id:
            raise ValueError(f'Dataset {dataset.dataset_id} is stored under key {dataset_id}')
        rd, _ = reduce_datablocks([[dataset]], reaction_prior)
        if len(rd) == 0:
            print(f'Dataset {dataset.dataset_id} skipped', file=sys.stderr)
            continue
        rd = rd[0][0]
        reduced_datasets[rd.dataset_id] = rd
    return reduced_datasets


# auxiliary functions for cli interface

def _read_input_file(filename: Optional[Path]=None) -> str:
    if filename is None:
        return None
    with open(filename, 'r') as f:
        return f.read()


def _read_data_file(filename: Optional[Path]=None) -> str:
    if filename is None:
        return None
    with open(filename, 'r') as f:
        return f.read()


def _require_input_cont(cont):
    if cont is None:
        print('Input file is missing, provide --input-file argument.')
        sys.exit(1)


def _require_data_cont(cont):
    if cont is None:
        print('Datasets file is missing, provide --data-file argument.')
        sys.exit(1)


if __name__ == '__main__':

    logging.basicConfig(
        level=logging.ERROR,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest='command', help='Available commands')

    extract_parser = subparsers.add_parser(
        'extract', help='Extract experimental data from GMDATA.CRD format'
    )
    extract_parser.add_argument(
        '--datasets', action='store_true', help='Extract the datasets from GMDATA.CRD'
    )
    extract_parser.add_argument(
        '--prior', action='store_true', help='Extract the reaction prior given in DAT.INP format'
    )
    extract_parser.add_argument(
        '--spectrum', action='store_true', help='Extract spectrum given in DAT.INP format'
    )
    extract_parser.add_argument(
        '--input-file', type=str, nargs='?', help='File with prior mesh and spectrum, usually named DAT.INP'
    )
    extract_parser.add_argument(
        '--data-file', type=str, nargs='?', help='File with experimental datasets, usually named GMDATA.CRD'
    )

    reduce_parser = subparsers.add_parser(
        'reduce', help='Reduce datasets to prior energy mesh'
    )
    reduce_parser.add_argument('filename', type=str, help='Input file with prior and datasets')

    args = parser.parse_args()

    if args.command == 'extract':
        input_cont = _read_input_file(args.input_file)
        data_cont = _read_data_file(args.data_file)
        combined_dict = {}
        if args.prior:
            _require_input_cont(input_cont)
            prior = extract_prior(input_cont)
            combined_dict['prior'] = prior
        if args.spectrum:
            _require_input_cont(input_cont)
            spectrum = extract_spectrum(input_cont)
            combined_dict['spectrum'] = spectrum
        if args.datasets:
            _require_data_cont(data_cont)
            datasets = extract_datasets(data_cont)
            combined_dict['datasets'] = datasets
        print(json.dumps(combined_dict, indent=2))
        sys.exit(0)

    elif args.command == 'reduce':
        filename = args.filename
        with open(filename, 'r') as f:
            stdin_cont = f.read()
        data_dict = json.loads(stdin_cont)
        prior = data_dict['prior']
        datasets = data_dict['datasets']
        reduce_result = reduce_datasets(datasets, prior)
        red_ds = {k: v.model_dump() for k, v in reduce_result.items()}
        print(json.dumps(red_ds, indent=2))
        sys.exit(0)
    else:
        raise ValueError('unknown command')
