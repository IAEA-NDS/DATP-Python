from typing import TextIO, Optional
import numpy as np
import math
from ...helpers import (
    fort_read,
    fort_write,
)
from ...constants import (
    ELIMINATION_BLOCK_INDICATION_STRING,
    DOWNWEIGHT_BLOCK_INDICATION_STRING,
    FISSION_SPECTRUM_BLOCK_INDICATION_STRING,
)
import logging


logger = logging.getLogger(__name__)


def _fort_write(file_handle, *args, **kwargs):
    if file_handle is None:
        return
    return fort_write(file_handle, *args, **kwargs)


def copy_gma_controls(
    prior_file_handle: TextIO,
    file_IO2: Optional[TextIO]=None,
    gma_file_handle: Optional[TextIO]=None
) -> dict:
    result = {}
    for K in range(10):
        format260 = '(A2,A2,A1,8I5)'
        KCO1, KCO2, KCO3, MC1, MC2, MC3, MC4, MC5, MC6, MC7, MC8 = \
            fort_read(prior_file_handle, format260, none_as=0.)

        logger.info(f'Reading control record {KCO1}')

        # exit loop if nothing more to read
        if KCO1.strip() == '':
            break

        _fort_write(gma_file_handle,  format260,
                   [KCO1, KCO2, KCO3, MC1, MC2, MC3, MC4, MC5, MC6, MC7, MC8])
        _fort_write(file_IO2,  format260,
                   [KCO1, KCO2, KCO3, MC1, MC2, MC3, MC4, MC5, MC6, MC7, MC8])
        if KCO1 == DOWNWEIGHT_BLOCK_INDICATION_STRING and MC2 == 10:
            # data set numbers selected for downweighting
            while True:
                format402 = '16I5'
                MSEP = np.empty((16,), dtype=int)
                MSEP[:] = fort_read(prior_file_handle, format402)
                _fort_write(file_IO2, format402, [MSEP])
                _fort_write(gma_file_handle, format402, [MSEP])
                if MSEP[0] == 0:
                    break

        elif KCO1 == FISSION_SPECTRUM_BLOCK_INDICATION_STRING and MC1 != 0:
            # fission spectrum
            energies = []
            spectrum_values = []
            while True:
                format404 = '(2E13.5)'
                AE, BS = fort_read(prior_file_handle, format404, none_as=0.)
                _fort_write(file_IO2, format404, [AE, BS])
                _fort_write(gma_file_handle, format404, [AE, BS])
                if AE == 0.0:
                    break
                energies.append(AE)
                spectrum_values.append(BS)

            result['spectrum'] = {
                'label': 'spectrum',
                'energies': energies,
                'spectrum_values': spectrum_values,
            }

        elif KCO1 == ELIMINATION_BLOCK_INDICATION_STRING:
            format408 = '(16i5)'
            format468 = "('Data Sets to be Excluded')"
            NEXL = fort_read(prior_file_handle, format408)
            _fort_write(gma_file_handle, format408, NEXL)
            _fort_write(file_IO2, format468, [None])
            _fort_write(file_IO2, format408, NEXL)
    return result
