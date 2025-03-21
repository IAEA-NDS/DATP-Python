import os
from .input_output import (
    copy_gma_controls,
    read_apriori,
    transfer_apriori_to_output_file,
    read_datablocks,
    GMADatabaseWriter
)
from .reduction import reduce_datablocks


def run_datp():

    basedir = '.'
    # OPEN(14,FILE='DAT.INP')
    prior_file_handle = open(os.path.join(basedir, 'DAT.INP'), 'r')
    # OPEN(15,FILE='DAT.LST')
    file_IO2 = open(os.path.join(basedir, 'DAT.LST'), 'w')
    # OPEN(12,FILE='GMDATA.CRD')
    expdata_file_handle = open(os.path.join(basedir, 'GMDATA.CRD'), 'r')
    # OPEN(13,FILE='DAT.RES')
    gma_file_handle = open(os.path.join(basedir, 'DAT.RES'), 'w')

    gmadb_writer = GMADatabaseWriter(gma_file_handle, file_IO2)

    copy_gma_controls(prior_file_handle, file_IO2, gma_file_handle)
    reaction_prior = read_apriori(prior_file_handle)
    transfer_apriori_to_output_file(file_IO2, gma_file_handle, reaction_prior)

    # read datablocks
    datablocks = read_datablocks(expdata_file_handle)
    reduced_datablocks = reduce_datablocks(datablocks, reaction_prior, file_IO2)

    # Produce the GMA Database file
    for datablock in reduced_datablocks:
        gmadb_writer.write_datablock_header()
        for dataset in datablock:
            gmadb_writer.write_dataset(dataset)
        gmadb_writer.write_datablock_trailer()
    gmadb_writer.write_trailer()

    prior_file_handle.close()
    file_IO2.close()
    expdata_file_handle.close()
    gma_file_handle.close()


if __name__ == '__main__':
    run_datp()
