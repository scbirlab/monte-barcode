"""Utilities used by monte-barcode."""

from collections.abc import Mapping
from functools import reduce
import operator
import os
import sys

import yaml

def _load_codons() -> Mapping:

    _data_path = os.path.join(os.path.dirname(__file__), 
                            'codons.yml')

    with open(_data_path, 'r') as f:
        codons = yaml.safe_load(f)

    codons['X'] = list(set(reduce(operator.add, list(codons.values()))))

    return codons


def _print_err(*args, **kwargs) -> None:

    return print(*args, **kwargs, file=sys.stderr)


def pprint_dict(x: Mapping, 
                message: str) -> None:
    
    key_val_str = (f'{key}: {val:.2f}' if isinstance(val, float) else f'{key}: {val}'
                   for key, val in x.items())

    _print_err(f'{message}:\n\t' + '\n\t'.join(key_val_str))
    
    return None

_CODONS = _load_codons()