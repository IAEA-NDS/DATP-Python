from .basemodels import (
    DatasetBase,
    ReactionPriorBase,
    FissionSpectrumBase,
)
import numpy as np


class ReactionPrior(ReactionPriorBase):
    pass


class Dataset(DatasetBase):

    def dict(self, *args, use_arrays=False, **kwargs):
        model_dict = super().dict(*args, **kwargs)
        if not use_arrays:
            return model_dict
        for key, value in model_dict.items():
            if isinstance(value, list) and key not in ('comments',):
                model_dict[key] = np.array(value)
        return model_dict


class FissionSpectrum(FissionSpectrumBase):
    pass
