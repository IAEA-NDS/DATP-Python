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


def simple_map(dest_key, source_key, cond):

    def forward_map(d):
        if cond == 'required':
            return {dest_key: d[source_key]}
        elif cond == 'optional':
            return {dest_key: d[source_key]} if source_key in d else {}
        else:
            raise ValueError('unsupported value of `cond`')

    def backward_map(d):
        if cond == 'required':
            return {source_key: d[dest_key]}
        elif cond == 'optional':
            return {source_key: d[dest_key]} if dest_key in d else {}
        else:
            raise ValueError('unsupported value of `cond`')

    return forward_map, backward_map


def identity_map(key, cond):
    """Identity map"""
    return simplemap(key, key, cond)


def oneway_func_map(key, func, cond):

    def forward_map(d):
        if cond == 'required':
            return {key: func(d)}
        elif cond == 'optional':
            try:
                return {key: func(d)}
            except Exception:
                return {}
    def backward_map(d):
        return {}

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
