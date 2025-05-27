import json
import os
from pathlib import Path
from copy import deepcopy
from .data_io.legacy.input_output import (
    copy_gma_controls,
)
from .data_io.legacy.data_input import (
    read_apriori,
    read_datablocks,
)
from .data_io.legacy.gma_output import (
    write_gmadb,
    write_prior,
)
from .reduction import reduce_datablocks
from .data_io.database import reduce_database
import argparse


def run_legacy_datp():

    basedir = '.'

    prior_file_handle = open(os.path.join(basedir, 'DAT.INP'), 'r')
    file_IO2 = open(os.path.join(basedir, 'DAT.LST'), 'w')
    expdata_file_handle = open(os.path.join(basedir, 'GMDATA.CRD'), 'r')
    gma_file_handle = open(os.path.join(basedir, 'DAT.RES'), 'w')

    copy_gma_controls(prior_file_handle, file_IO2, gma_file_handle)
    reaction_prior = read_apriori(prior_file_handle)
    write_prior(file_IO2, gma_file_handle, reaction_prior)

    datablocks = read_datablocks(expdata_file_handle)
    reduced_datablocks, auxinfo_blocks = (
        reduce_datablocks(datablocks, reaction_prior, file_IO2)
    )
    write_gmadb(gma_file_handle, file_IO2, reduced_datablocks, auxinfo_blocks)

    prior_file_handle.close()
    file_IO2.close()
    expdata_file_handle.close()
    gma_file_handle.close()


def _reduce_database(gmadb):
    reaction_prior = map_priorblocks(gmadb['prior'], direction='forward')
    datablocks = map_datablocks(gmadb['datablocks'], direction='forward')
    reduced_datablocks = []
    # NOTE: Fission spectrum is assumed to come last in prior.
    #       If not true, prior indexing would produce garbage.
    assert reaction_prior[-1]['type'] == 'legacy-fission-spectrum'

    new_datablocks = []
    for datablock in datablocks:
        new_datablock = deepcopy(datablock)
        reduced_datablock = _ = (
            reduce_datablocks([datablock['datasets']], reaction_prior)
        )
        if len(reduced_datablock) == 1:
            new_datablock['datasets'] = reduced_datablock[1]
            new_datablocks.append(new_datablock)

    new_gmadb = deepcopy(gmadb)
    new_gmadb['datablocks'] = new_datablocks
    return new_gmadb


def run_datp(dbfile_in, dbfile_out):
    """Read GMA database file and output reduced one."""
    with open(dbfile_in, 'r') as f:
        gmadb = json.load(f)
    new_gmadb = _reduce_database(gmadb)
    with open(dbfile_out, 'w') as f:
        json.dump(new_gmadb, dbfile_out, indent=2)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--legacy', action='store_true', help='datpy behaves as if it were the DATP Fortran code')
    parser.add_argument('--input', help='the database with datasets to be reduced')
    parser.add_argument('--output', help='reduced database is written to this file')
    args = parser.parse_args()

    if args.legacy:
        if args.input or args.output:
            parser.error('--legacy cannot be used with --input or --output')
    else:
        if not args.input or not args.output:
            parser.error('--input and --output are required unless --legacy is specified')

    if args.legacy is True:
        run_legacy_datpy()
    else:
        input_file = Path(args.input)
        output_file = Path(args.output)
        run_datp(input_file, output_file)
