import os
import io
import sys
import json
from typing import Dict
from .input_output import (
    copy_gma_controls,
)
from .data_input import (
    read_apriori,
    read_datablocks,
)
from .input_output import (
    copy_gma_controls,
)
from .gma_output import (
    write_gmadb,
    write_prior,
)
from .reduction import reduce_datablocks
import argparse
import logging
from .datamodels.models import Dataset
from .reduction import reduce_datablocks


logger = logging.getLogger(__name__)


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


def reduce_datasets(datasets: Dict[int, dict]) -> dict:
    reduced_datasets = {}
    for dataset_id, raw_dataset in datasets.items():
        dataset = Dataset(**raw_dataset)
        if int(dataset_id) != dataset.dataset_id:
            raise ValueError(f'Dataset {dataset.dataset_id} is stored under key {dataset_id}')
        rd, _ = reduce_datablocks([[dataset]])
        reduced_datasets[rd.dataset_id] = rd.dict()
    return reduced_datasets


if __name__ == '__main__':

    logging.basicConfig(
        level=logging.ERROR,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

    stdin_cont = sys.stdin.read()

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

    args = parser.parse_args()

    if args.command == 'extract':
        combined_dict = {}
        if args.prior:
            prior = extract_prior(stdin_cont)
            combined_dict['prior'] = prior
        if args.spectrum:
            spectrum = extract_spectrum(stdin_cont)
            combined_dict['spectrum'] = spectrum
        if args.datasets:
            datasets = extract_datasets(stdin_cont)
            combined_dict['datasets'] = datasets
        print(json.dumps(combined_dict, indent=2))
        sys.exit(0)
    elif args.command == 'reduce':
        datasets = json.loads(stdin_cont)
        print(json.dumps(reduce_datasets(datasets), indent=2))
        sys.exit(0)
    else:
        raise ValueError('unknown command')
