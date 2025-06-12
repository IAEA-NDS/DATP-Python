def is_empty(arr):
    if isinstance(arr, list):
        return all(is_empty(a) for a in arr)
    return arr is None


def transpose2d(arr):
    lens = [len(a) for a in arr]
    if not all(lens[0] == l for l in lens):
        raise ValueError('lists must be of same length')
    if not all([isinstance(a1, (float, int)) for a1 in a2] for a2 in arr):
        raise TypeError('inner lists must only contain numbers')
    lena = lens[0]
    return [[a[i] for a in arr] for i in range(lena)]


def apply_transpose2d(arr, transpose):
    if transpose:
        return transpose2d(arr)
    return arr


def simple_augmentation(key, value, direction):

    def forward_map(d):
        if direction == 'forward':
            return {key: value}
        else:
            return {}
        
    def backward_map(d):
        if direction == 'backward':
            return {key: value}
        else:
            return {}

    return forward_map, backward_map


def simple_map(dest_key, source_key, cond, transpose=False, drop_empty=True):

    def forward_map(d):
        if cond == 'required':
            val = d[source_key]
            if is_empty(val) and drop_empty:
                return {}
            return {dest_key: apply_transpose2d(val, transpose)}
        elif cond == 'optional':
            if source_key not in d:
                return {}
            val = d[source_key]
            if is_empty(val) and drop_empty:
                return {}
            return {dest_key: apply_transpose2d(val, transpose)}
        else:
            raise ValueError('unsupported value of `cond`')

    def backward_map(d):
        if cond == 'required':
            val = d[dest_key]
            if is_empty(val) and drop_empty:
                return {}
            return {source_key: apply_transpose2d(val, transpose)}
        elif cond == 'optional':
            if dest_key not in d:
                return {}
            val = d[dest_key]
            if is_empty(val) and drop_empty:
                return {}
            return {source_key: apply_transpose2d(d[dest_key], transpose)}
        else:
            raise ValueError('unsupported value of `cond`')

    return forward_map, backward_map


def identity_map(key, cond, transpose=False, drop_empty=True):
    """Identity map"""
    return simple_map(key, key, cond, transpose, drop_empty)


def oneway_func_map(key, func, cond, direction='forward'):

    def custom_map(d):
        if cond == 'required':
            return {key: func(d)}
        elif cond == 'optional':
            try:
                return {key: func(d)}
            except Exception:
                return {}

    if direction == 'forward':
        forward_map = custom_map
        backward_map = lambda d: {}
    else:
        forward_map = lambda d: {}
        backward_map = custom_map

    return forward_map, backward_map


def update_dict(dest_dict, source_dict, mapping, direction="forward"):
    if direction == "forward":
        dest_dict.update(mapping[0](source_dict))
    elif direction == "backward":
        dest_dict.update(mapping[1](source_dict))
    else:
        raise ValueError('unsupported value of `direction`')


def map_dict(source_dict: dict, mappings: list, direction='forward'):
    dest_dict = {}
    for curmap in mappings:
        update_dict(dest_dict, source_dict, curmap, direction)
    return dest_dict
