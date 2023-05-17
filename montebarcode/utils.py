"""Utilities used by monte-barcode."""

from collections.abc import Mapping
from functools import reduce
import operator
import os
import sys

import yaml

_data_path = os.path.join(os.path.dirname(__file__), 
                          'codons.yml')

with open(_data_path, 'r') as f:
    _CODONS = yaml.safe_load(f)

_CODONS['X'] = list(set(reduce(operator.add, list(_CODONS.values()))))

def pprint_dict(x: Mapping, 
                message: str) -> None:
    
    key_val_str = (f'{key}: {val:.2f}' if isinstance(val, float) else f'{key}: {val}'
                   for key, val in x.items())

    print(f'{message}:\n\t' + '\n\t'.join(key_val_str),
          file=sys.stderr)
    
    return None