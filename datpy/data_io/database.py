from .database_mappings import (
    map_priorblocks,
    map_datablocks,
)


def assemble_database(
    reaction_prior: ReactionPrior, datablocks: list[Dataset], fission_spectrum: FissionSpectrum
) -> dict:
    reaction_prior = reaction_prior.dict()
    datablocks = [blck.dict() for blck in datablocks]
    fission_spectrum = fission_spectrum.dict()

    new_db = {}
    new_db['prior'] = map_priorblocks(reaction_prior, direction='backward')
    new_db['prior'].append(map_priorblock(fission_spectrum, direction='backward'))

    mapped_datablocks = map_datablocks(datablocks, direction='backward')
    # add the outer shell to the datablocks
    new_datablocks = []
    for datablock in mapped_datablocks:
        new_datablock = {
            'type': 'legacy-experiment-datablock', 
            'datasets': datablock
        }
        if any('ECOR' in ds for ds in datablock[:-1]):
            raise ValueError(
                'ECOR field only allowed in last dataset of datablock'
            )
        if 'ECOR' in datablock[-1]: 
            new_datablock['ECOR'] = datablock[-1].pop('ECOR') 
        new_datablocks.append(new_datablock)

    new_db['datablocks'] = new_datablocks
    return new_db


def convert_database_to_internal_format(database):
    reaction_prior = map_priorblocks(gmadb['prior'], direction='forward')
    datablocks = map_datablocks(gmadb['datablocks'], direction='forward')
    # extract the fissions spectrum from prior
    spectrum_idx = [
        i for i, d in enumerate(reaction_prior) if d['type'] == 'legacy-fission-spectrum'
    ]
    if len(spectrum_idx) > 1:
        raise ValueError('Only one fission spectrum allowed')
    fission_spectrum = None
    if len(spectrum_idx) == 1:
        fission_spectrum = reaction_prior[spectrum_idx[0]]

    return (
        reaction_prior,
        datablocks,
        fission_spectrum,
    )
