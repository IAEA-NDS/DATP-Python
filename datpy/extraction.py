import os
import io
import sys
import json
from .input_output import (
    copy_gma_controls,
)
from .data_input import (
    read_apriori,
    read_datablocks,
)
from .gma_output import (
    write_gmadb,
    write_prior,
)
from .reduction import reduce_datablocks
import argparse
import logging


logger = logging.getLogger(__name__)


def extract_datasets(gmdata_crd_content: str):

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



if __name__ == '__main__':

    logging.basicConfig(
        level=logging.ERROR,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

    parser = argparse.ArgumentParser()
    parser.add_argument('filename')

    args = parser.parse_args()
    fn = args.filename
    with open(fn, 'r') as f:
        cont = f.read()

    datasets = extract_datasets(cont)
    print(json.dumps(datasets, indent=2))
