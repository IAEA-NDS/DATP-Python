from .constants import (
    MAXF,
    FORMAT200,
    FORMAT250,
    FORMAT290,
)
from .helpers import (
    fort_write
)


def write_gmadb_file(gma_file_handle, file_IO2, reduced_datablocks):
    for datablock in reduced_datablocks:
        write_datablock_header(gma_file_handle, file_IO2)
        for dataset in datablock:
            write_dataset(gma_file_handle, file_IO2, dataset)
        write_datablock_trailer(gma_file_handle, file_IO2)
    write_file_trailer(gma_file_handle, file_IO2)


def write_datablock_header(gma_file_handle, file_IO2):
    format251 = '(4HBLCK,1X,2I5)'
    fort_write(gma_file_handle, format251, [0, 0])
    fort_write(file_IO2, format251, [0, 0])


def write_datablock_trailer(gma_file_handle, file_IO2):
    fort_write(gma_file_handle, FORMAT250, [0, 0])
    fort_write(file_IO2, FORMAT250, [0, 0])


def write_file_trailer(gma_file_handle, file_IO2): 
    format256 = '(4HEND*,1X,2I5)'
    fort_write(gma_file_handle, format256, [0, 0])
    fort_write(file_IO2, format256, [0, 0])


def write_dataset(gma_file_handle, file_IO2, dataset):
    ds = dataset
    # find nr of CS involved
    for MN in range(1, 5):  # 66
        pyMN = MN - 1
        if ds.reaction_ids[pyMN] == 0:
            break
    NNN = MN - 1

    format253 = '(5HDATA ,9I5)'
    fort_write(gma_file_handle, format253, [ds.dataset_id, ds.quantity_type, ds.cormat_dim, NNN, ds.reaction_ids[0:4]])
    fort_write(file_IO2, format253, [ds.dataset_id, ds.quantity_type, ds.cormat_dim, NNN, ds.reaction_ids[0:4]])

    format254 = '(3I5,A28,8X,A20)'
    fort_write(gma_file_handle, format254, [ds.year, ds.NQ, ds.NCST_size, ds.author, ds.pubref])
    fort_write(file_IO2, format254, [ds.year, ds.NQ, ds.NCST_size, ds.author, ds.pubref])

    if ds.quantity_type not in (2, 4, 8, 9):
        # normalization uncertainties
        format261 = '(10F5.1,10I3)'
        fort_write(gma_file_handle, format261, [ds.ENF[0:10], ds.NENF[0:10]])
        fort_write(file_IO2, format261, [ds.ENF[0:10], ds.NENF[0:10]])

    # energy dep. unc. parameters
    format262 = '(3F5.2,I3)'
    # NOTE: this loop is implicit in the write statement
    #       of the fortran code
    for K in range(11):
        fort_write(file_IO2, format262, [ds.EPA[0:3, K], ds.NETG[K]])
        fort_write(gma_file_handle, format262, [ds.EPA[0:3, K], ds.NETG[K]])
    if ds.NCST_size != 0:
        # cross correlations
        format263 = '(I5,20I3)'
        format293 = '(10F5.1)'
        for K in range(ds.NCST_size):  # 83
            # NOTE: during flattening in fort_write first index should
            #       change fastest
            fort_write(gma_file_handle, format263, [ds.NCST[K], ds.NEC[:, :, K]])
            fort_write(file_IO2, format263, [ds.NCST[K], ds.NEC[:, :, K]])
            fort_write(gma_file_handle, format293, [ds.FCFC[0:10, K]])
            fort_write(file_IO2, format293, [ds.FCFC[0:10, K]])

    if ds.quantity_type == 6:
        # fission spectrum average data set
        fort_write(file_IO2, FORMAT200, [ds.energies[0], ds.measured_values[0], ds.uncertainties[0:12, 0]])
        fort_write(gma_file_handle, FORMAT200, [ds.energies[0], ds.measured_values[0], ds.uncertainties[0:12, 0]])
        fort_write(file_IO2, FORMAT200, [0, 0, ds.uncertainties[0:12, MAXF-1]])
        fort_write(gma_file_handle, FORMAT200, [0, 0, ds.uncertainties[0:12, MAXF-1]])
        return

    # write out more for data types different from NT==6
    format5173 = "(/' ENERGY/MEV  VALUE       UNCERTAINTIES                     RATIO TO APRIORI'/)"
    fort_write(file_IO2, format5173, [None])
    for k in range(ds.NO):
        EEE = ds.energies[k]
        QQQ = ds.measured_values[k]
        DIF = ds.DIF[k]
        fort_write(gma_file_handle, FORMAT200, [EEE, QQQ, ds.uncertainties[0:12, k]])
        fort_write(file_IO2, FORMAT290, [EEE, QQQ, ds.uncertainties[0:12, k], DIF])

    # end of data set
    fort_write(gma_file_handle, FORMAT200, [0, 0, ds.uncertainties[0:12, MAXF-1]])
    fort_write(file_IO2, FORMAT200, [0, 0, ds.uncertainties[0:12, MAXF-1]])

    if ds.cormat_dim == 0:
        return

    format6114 = '(1X,10F8.5)'
    format6115 = '(10F8.5)'
    for KL in range(ds.cormat_dim):  # 6113
        fort_write(file_IO2, format6114, [ds.cormat[KL, :(KL+1)]])
        fort_write(gma_file_handle, format6115, [ds.cormat[KL, :(KL+1)]])
