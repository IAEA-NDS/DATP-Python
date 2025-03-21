import os
from .input_output import (
    copy_gma_controls,
    read_apriori,
    transfer_apriori_to_output_file,
    read_datablocks,
)
from .gma_output import write_gmadb_file
from .reduction import reduce_datablocks


def run_datp():

    basedir = '.'

    prior_file_handle = open(os.path.join(basedir, 'DAT.INP'), 'r')
    file_IO2 = open(os.path.join(basedir, 'DAT.LST'), 'w')
    expdata_file_handle = open(os.path.join(basedir, 'GMDATA.CRD'), 'r')
    gma_file_handle = open(os.path.join(basedir, 'DAT.RES'), 'w')

    copy_gma_controls(prior_file_handle, file_IO2, gma_file_handle)
    reaction_prior = read_apriori(prior_file_handle)
    transfer_apriori_to_output_file(file_IO2, gma_file_handle, reaction_prior)

    datablocks = read_datablocks(expdata_file_handle)
    reduced_datablocks = reduce_datablocks(datablocks, reaction_prior, file_IO2)
    write_gmadb_file(gma_file_handle, file_IO2, reduced_datablocks)

    prior_file_handle.close()
    file_IO2.close()
    expdata_file_handle.close()
    gma_file_handle.close()


if __name__ == '__main__':
    run_datp()
