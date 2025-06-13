import json
import os
from typing import Optional
from pathlib import Path
from copy import deepcopy
from .datamodels.models import (
    Dataset,
    ReactionPrior,
)
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
from .data_io.database_mappings import (
    map_priorblocks,
    map_priorblock,
    map_datablocks,
)
from .reduction import reduce_datablocks
import argparse


def run_legacy_datp(dbfile_out: Optional[str]=None, do_reduce=True):
    """Run datpy in legacy DATP mode

    Parameters
    ----------
    dbfile_out : str
        Write resulting data to a JSON file.
        If ``None``, store the results in the same
        file and using the same format as legacy DATP.
    do_reduce : bool
        If ``True``, perform the data reduction, otherwise
        just transfer the data as they are to the output files.

    Returns
    -------
    None
        This function returns nothing.
    """
    basedir = '.'

    if dbfile_out is None:
        file_IO2 = open(os.path.join(basedir, 'DAT.LST'), 'w')
        gma_file_handle = open(os.path.join(basedir, 'DAT.RES'), 'w')
    else:
        file_IO2 = None
        gma_file_handle = None

    prior_file_handle = open(os.path.join(basedir, 'DAT.INP'), 'r')
    expdata_file_handle = open(os.path.join(basedir, 'GMDATA.CRD'), 'r')

    # The `copy_gma_controls` function does not write anything
    # if None is provided as output file handles. However, the function
    # still needs to be called as it advances the file pointer
    # of the input file containing the reaction prior.
    spectrum_dict = copy_gma_controls(prior_file_handle, file_IO2, gma_file_handle)
    reaction_prior = read_apriori(prior_file_handle)
    datablocks = read_datablocks(expdata_file_handle)

    if do_reduce:
        reduced_datablocks, auxinfo_blocks = (
            reduce_datablocks(datablocks, reaction_prior, file_IO2)
        )
    else:
        reduced_datablocks = datablocks

    if dbfile_out is None:
        if not do_reduce:
            # NOTE: The reason is the calculation of the `auxinfo_blocks`
            #       data structure during reduction, which is also written
            #       to the legacy output file.
            raise NotImplementedError(
                "Writing the results to output files in legacy format "
                "is only supported in combination with data reduction."
            )
        write_prior(file_IO2, gma_file_handle, reaction_prior)
        write_gmadb(gma_file_handle, file_IO2, reduced_datablocks, auxinfo_blocks)
        file_IO2.close()
        gma_file_handle.close()
    else:
        # map prior and fission spectrum to gmapy format
        reaction_prior_out = map_priorblocks(
            reaction_prior.model_dump(), direction='backward', do_reduce=do_reduce)
        reaction_prior_out.append(
            map_priorblock(spectrum_dict['spectrum'], direction='backward', do_reduce=do_reduce)
        )
        # map datablocks to gmapy format
        red_db = [[ds.model_dump() for ds in db] for db in reduced_datablocks]
        datablocks_out = map_datablocks(red_db, direction='backward')
        # put everything together and write to json file
        dbout = {'prior': reaction_prior_out, 'datablocks': datablocks_out}
        with open(dbfile_out, 'w') as f:
            json.dump(dbout, f, indent=2)

    prior_file_handle.close()
    expdata_file_handle.close()


def _reduce_database(gmadb):
    if gmadb['prior'][-1]['type'] != 'legacy-fission-spectrum':
        raise ValueError(
            'The last list item in `prior` is expected to be of type `legacy-fission-spectrum`'
        )
    spectrum_raw = gmadb['prior'][-1]
    reaction_prior_raw = gmadb['prior'][:-1]

    reaction_prior = map_priorblocks(
        reaction_prior_raw, direction='forward', do_reduce=True
    )
    reaction_prior_internal = ReactionPrior(reaction_prior)

    datablocks = map_datablocks(gmadb['datablocks'], direction='forward')
    datablocks_internal = [[Dataset(**ds) for ds in b] for b in datablocks]

    new_datablocks_out = []
    for datablock in datablocks_internal:
        reduced_datablock, _ = (
            reduce_datablocks([datablock], reaction_prior_internal)
        )
        if len(reduced_datablock) == 1:
            new_datablock = [ds.model_dump() for ds in reduced_datablock[0]]
            new_datablock_out = map_datablocks([new_datablock], direction='backward')[0]
            new_datablocks_out.append(new_datablock_out)

    # NOTE: The prior mesh in the input file to DATP (DAT.INP or JSON file) is
    #       NOT identical to the one in the output file. The first
    #       and last mesh point are removed in the output file.
    #       The `reduce_datablocks` routine expects the full prior mesh
    #       and the prior cross section dictionary with field names
    #       as used internally in `datpy` described in datamodels/schemas.py.
    #       Therefore, afterwards, we need to convert back and the `do_reduce` flag
    #       takes care of removing the point.
    reaction_prior_out = map_priorblocks(
        reaction_prior, direction='backward', do_reduce=True
    )

    reaction_prior_out = map_priorblocks(
        reaction_prior_internal.model_dump(), direction='backward', do_reduce=True
    )
    new_gmadb = {
        'prior': reaction_prior_out + [spectrum_raw],
        'datablocks': new_datablocks_out
    }
    return new_gmadb


def run_datp(dbfile_in, dbfile_out):
    """Read GMA database file and output reduced one."""
    with open(dbfile_in, 'r') as f:
        gmadb = json.load(f)

    new_gmadb = _reduce_database(gmadb)
    with open(dbfile_out, 'w') as f:
        json.dump(new_gmadb, f, indent=2)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--legacy', action='store_true', help='datpy behaves as if it were the DATP Fortran code')
    parser.add_argument('--no-reduce', action='store_true', help='do not perform reduction')
    parser.add_argument('--input', help='the database with datasets to be reduced')
    parser.add_argument('--output', help='reduced database is written to this file')
    args = parser.parse_args()

    if args.legacy:
        if args.input:
            parser.error('--legacy cannot be used with --input')
    else:
        if args.no_reduce:
            parser.error('The --no-reduce argument can only be used together with --legacy')
        if not args.input or not args.output:
            parser.error('--input and --output are required unless --legacy is specified')

    if args.legacy is True:
        run_legacy_datp(dbfile_out=args.output, do_reduce=(not args.no_reduce))
    else:
        input_file = Path(args.input)
        output_file = Path(args.output)
        run_datp(input_file, output_file)
