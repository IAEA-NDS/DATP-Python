from .basemodels import DatasetBase
import numpy as np


class Dataset(DatasetBase):

    def dict(self, *args, use_arrays=False, **kwargs):
        model_dict = super().dict(*args, **kwargs)
        for key, value in model_dict.items():
            if isinstance(value, list) and key not in ('comments',):
                model_dict[key] = np.array(value)

        return model_dict
