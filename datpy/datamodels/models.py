from .basemodels import (
    DatasetBase,
    ReactionPriorBase,
    FissionSpectrumBase,
)
import numpy as np


class ReactionPrior(ReactionPriorBase):
    pass


class Dataset(DatasetBase):

    def _model_dump(self, model_dict, use_arrays):
        if not use_arrays:
            return model_dict
        for key, value in model_dict.items():
            if isinstance(value, list) and key not in ('comments',):
                model_dict[key] = np.array(value)
        return model_dict

    def dict(self, *args, use_arrays=False, **kwargs):
        model_dict = super().dict(*args, **kwargs)
        return self._model_dump(model_dict, use_arrays)

    def model_dump(self, *args, use_arrays=False, **kwargs):
        model_dict = super().model_dump(*args, **kwargs)
        return self._model_dump(model_dict, use_arrays)


class FissionSpectrum(FissionSpectrumBase):
    pass
