import numpy as np
from ...constants import (
    FORMAT200,
    FORMAT250,
    FORMAT290,
    SHOULD_TEST_OUTPUT,
)
from ...helpers import (
    Bunch,
    fort_write,
)


def write_gmadb(gma_file_handle, file_IO2, reduced_datablocks, auxinfo_blocks):
    for datablock, auxinfo in zip(reduced_datablocks, auxinfo_blocks):
        write_datablock_header(gma_file_handle, file_IO2)
        for dataset in datablock:
            write_dataset(gma_file_handle, file_IO2, dataset, auxinfo)
        write_datablock_trailer(gma_file_handle, file_IO2)
    write_file_trailer(gma_file_handle, file_IO2)


def write_prior(file_IO2, gma_file_handle, reaction_prior):
    reaction_prior = reaction_prior.dict()
    prior_label =  [p['label'] for p in reaction_prior.values()]
    prior_energy_mesh = [p['energies'] for p in reaction_prior.values()]
    prior_cross_section = [p['cross_sections'] for p in  reaction_prior.values()]

    num_reactions = len(prior_energy_mesh)
    ITOT = 0
    for K in range(num_reactions):
        ITOT = ITOT + len(prior_energy_mesh[K])
    ITOT = ITOT - 2*num_reactions

    format264 = '(5HAPRI ,2I5)'
    fort_write(gma_file_handle, format264, [ITOT, num_reactions])
    fort_write(file_IO2, format264, [ITOT, num_reactions])

    format3731 = "('cross section',i5,' number',i5,A16)"
    format101 = '(2X,I4,2X,2E12.4)'
    format109 = '(8x,2e10.4)'
    format100 = '(2E14.6)' if SHOULD_TEST_OUTPUT else '(2E10.4)'
    format99 = '(A16)'  # original: (8A2)
    for L in range(num_reactions):
        NOR = len(prior_energy_mesh[L]) - 1
        NOR2 = NOR - 1
        fort_write(gma_file_handle, format99, [prior_label[L]])
        fort_write(file_IO2, format99, [prior_label[L]])
        fort_write(None, format3731, [L+1, NOR2, prior_label[L]])

        for K in range(1, NOR):
            fort_write(file_IO2, format101, [K, prior_energy_mesh[L][K], prior_cross_section[L][K]])
            fort_write(gma_file_handle, format100, [prior_energy_mesh[L][K], prior_cross_section[L][K]])
        fort_write(file_IO2, format109, [0, 0])
        fort_write(gma_file_handle, format100, [0, 0])


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


def write_dataset(gma_file_handle, file_IO2, dataset, auxinfo):
    ds = Bunch(dataset.dict(use_arrays=True))
    # find number of CS involved
    NNN = ds.num_reaction_ids

    if NNN > 3:
        raise ValueError('All quantity types have not more than three reaction ids')
    format253 = '(5HDATA ,8I5)'
    fort_write(gma_file_handle, format253, [ds.dataset_id, ds.quantity_type, ds.cormat.shape[0], NNN, ds.reaction_ids[0:3], ds.NNCOX])
    fort_write(file_IO2, format253, [ds.dataset_id, ds.quantity_type, ds.cormat.shape[0], NNN, ds.reaction_ids[0:3], ds.NNCOX])

    format254 = '(3I5,A28,8X,A20)'
    fort_write(gma_file_handle, format254, [ds.year, ds.tag, len(ds.NCST), ds.author, ds.pubref])
    fort_write(file_IO2, format254, [ds.year, ds.tag, len(ds.NCST), ds.author, ds.pubref])

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
    if len(ds.NCST) != 0:
        # cross correlations
        format263 = '(I5,20I3)'
        format293 = '(10F5.1)'
        for K in range(len(ds.NCST)):  # 83
            # NOTE: during flattening in fort_write first index should
            #       change fastest
            fort_write(gma_file_handle, format263, [ds.NCST[K], np.transpose(ds.NEC[:, :, K])])
            fort_write(file_IO2, format263, [ds.NCST[K], ds.NEC[:, :, K]])
            fort_write(gma_file_handle, format293, [ds.FCFC[0:10, K]])
            fort_write(file_IO2, format293, [ds.FCFC[0:10, K]])

    if ds.quantity_type == 6:
        # fission spectrum average data set
        fort_write(file_IO2, FORMAT200, [ds.energies[0], ds.measured_values[0], ds.uncertainties[0:12, 0]])
        fort_write(gma_file_handle, FORMAT200, [ds.energies[0], ds.measured_values[0], ds.uncertainties[0:12, 0]])
        reduced_uncertainties = [0.0]*12
        fort_write(file_IO2, FORMAT200, [0, 0, reduced_uncertainties])
        fort_write(gma_file_handle, FORMAT200, [0, 0, reduced_uncertainties])
        return

    # write out more for data types different from NT==6
    format5173 = "(/' ENERGY/MEV  VALUE       UNCERTAINTIES                     RATIO TO APRIORI'/)"
    fort_write(file_IO2, format5173, [None])
    for k in range(len(ds.measured_values)):
        EEE = ds.energies[k]
        QQQ = ds.measured_values[k]
        DIF = auxinfo[ds.dataset_id]['DIF'][k]
        fort_write(gma_file_handle, FORMAT200, [EEE, QQQ, ds.uncertainties[0:12, k]])
        fort_write(file_IO2, FORMAT290, [EEE, QQQ, ds.uncertainties[0:12, k], DIF])

    # end of data set
    reduced_uncertainties = auxinfo[ds.dataset_id]['reduced_uncertainties']
    fort_write(gma_file_handle, FORMAT200, [0, 0, reduced_uncertainties])
    fort_write(file_IO2, FORMAT200, [0, 0, reduced_uncertainties])

    if not hasattr(ds, 'cormat') or ds.cormat.shape[0] == 0:
        return

    format6114 = '(1X,10F8.5)'
    format6115 = '(10F8.5)'
    for KL in range(ds.cormat.shape[0]):  # 6113
        fort_write(file_IO2, format6114, [ds.cormat[KL, :(KL+1)]])
        fort_write(gma_file_handle, format6115, [ds.cormat[KL, :(KL+1)]])
