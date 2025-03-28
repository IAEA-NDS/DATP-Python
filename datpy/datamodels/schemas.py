__all__ = ['schema_list']

schema_list = []

# auxiliary info

_intstr_regex = '^[0-9]+$'

# recurring properies

def _int_property(min_value=None, max_value=None):
    ret = {'type': 'integer'}
    if min_value is not None:
        ret['minimum'] = min_value
    if max_value is not None:
        ret['maximum'] = max_value
    return ret


def _float_property(min_value=None, max_value=None):
    ret = {'type': 'number'}
    if min_value is not None:
        ret['minimum'] = min_value
    if max_value is not None:
        ret['maximum'] = max_value
    return ret


def _str_property(min_length=None, max_length=None):
    ret = {'type': 'string'}
    if min_length is not None:
        ret['minLength'] = min_length
    if max_length is not None:
        ret['maxLength'] = max_length
    return ret


def _enum_property(choice, dtype):
    return {
        'type': dtype,
        'enum': choice,
    }


def _array_property(elem_property=None, min_items=None, max_items=None):
    ret = {'type': 'array'}
    if elem_property is not None:
        ret['items'] = elem_property
    if min_items is not None:
        ret['minItems'] = min_items
    if max_items is not None:
        ret['maxItems'] = max_items
    return ret


def _multidim_array_property(elem_property, sizes):
    elem_property = _array_property(elem_property, sizes[-1], sizes[-1])
    if len(sizes) > 1:
        return _multidim_array_property(elem_property, sizes[:-1])
    return elem_property


reaction_prior = {
    '$schema': 'https://json-schema.org/draft/2020-12/schema',
    'version': '0.0.1',
    'title': 'ReactionPriorBase',
    'description': 'List of energy-dependent cross sections with prior values',
    'type': 'object',
    'patternProperties': {
        _intstr_regex: {
            'type': 'object',
            'properties': {
                'label': _str_property(16, 16),
                'reaction_id': _int_property(1, None),
                'energies': _array_property(_float_property()),
                'cross_sections': _array_property(_float_property()),
            },
            'required': ['label', 'energies', 'cross_sections']
        }
    },
    'additionalProperties': False,
}
schema_list.append(reaction_prior)


dataset_schema = {
    '$schema': 'https://json-schema.org/draft/2020-12/schema',
    'version': '0.0.1',
    'title': 'DatasetBase',
    'description': 'Dataset with experimental measurements and uncertainties',
    'type': 'object',
    'properties': {
        'dataset_id': _int_property(0, 9999),
        'year': _int_property(1900, None),
        'author': _str_property(None, 28),
        'pubref': _str_property(None, 20),
        'tag': _int_property(0, 99),
        'quantity_type': _int_property(0, 9),
        'comments': _array_property(_str_property(None, 80)),
        'num_reaction_ids': _int_property(1, 5),
        'reaction_ids': _array_property(_int_property(0, None), 1, 5),
        'ENF': _array_property(_float_property(), 10, 10),
        'NENF': _array_property(_int_property(0, 999), 10, 10),
        'EPA': _multidim_array_property(_float_property(), [3, 11]),
        'NETG': _array_property(_int_property(0, 999), 11, 11),
        'energies': _array_property(_float_property()),
        'measured_values': _array_property(_float_property()),
        'uncertainties': _multidim_array_property(_float_property(), [12, None]),
        'NCST': _array_property(_int_property(0, 9999)),
        'NEC': _multidim_array_property(_int_property(0, 21), [2, 10, None]),
        'FCFC': _multidim_array_property(_float_property(), [10, None]),
        'cormat': _multidim_array_property(_float_property(-1.0, 1.0), [None, None]),
    },
    'required': [
       'dataset_id', 'year', 'author', 'pubref', 'tag', 'quantity_type',
       'comments', 'num_reaction_ids', 'reaction_ids', # not required: ENF, NENF
       'EPA', 'NETG', 'energies', 'measured_values', 'uncertainties',
       'NCST', 'NEC', 'FCFC',  # not required: cormat
    ]
}
schema_list.append(dataset_schema)
