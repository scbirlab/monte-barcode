"""Functions for generating random barcodes."""

from __future__ import annotations

from typing import Union
from collections.abc import Generator, Iterable, Mapping, Sequence
from collections import Counter
from functools import reduce
from itertools import groupby, product
import operator
from random import choices, sample

import streq as sq

from .utils import _CODONS

def _tmat(x: Iterable[str]) -> Mapping:
    
    xy = sorted(zip(*x), 
                key=lambda x: x[0])
    
    counter = {key: Counter(following for _, following in group) 
               for key, group in groupby(xy, lambda x: x[0])}

    return [{key: (tuple(c.keys()), tuple(c.values())) for key, c in counter.items()}]


def transition_matrix(x: Sequence[str]) -> Sequence[Mapping]:

    """Generate transition frequencies from one item to the next in a sequence.
    
    Counts the occurence of the next letter conditioned on the preceding letter.

    Parameters
    ----------
    x : Sequence[str]
        List of strings to take transition frequencies from.

    Returns
    -------
    tuple
        A length-n tuple, where n is the minimum length of x. Each item is a 
        2-tuple containing the next possible letters and their frequencies.

    Examples
    --------
    >>> transition_matrix(['ATC', 'ATG'])  # doctest: +NORMALIZE_WHITESPACE
    ({None: (('A',), (2,))}, {'A': (('T',), (2,))}, {'T': (('C', 'G'), (1, 1))})
    >>> transition_matrix(['ATC', 'CTG'])  # doctest: +NORMALIZE_WHITESPACE
    ({None: (('A', 'C'), (1, 1))}, {'A': (('T',), (1,)), 'C': (('T',), (1,))}, {'T': (('C', 'G'), (1, 1))})
    >>> transition_matrix(['ATC', 'CAG'])  # doctest: +NORMALIZE_WHITESPACE
    ({None: (('A', 'C'), (1, 1))}, {'A': (('T',), (1,)), 'C': (('A',), (1,))}, {'A': (('G',), (1,)), 'T': (('C',), (1,))})
    
    """

    initial = Counter(_x[0] for _x in x)
    initial = [{None: (tuple(initial.keys()), tuple(initial.values()))}]

    preceding = zip(*x)
    following = zip(*(_x[1:] for _x in x))

    return tuple(reduce(operator.add, map(_tmat, zip(preceding, following)), initial))


def codon_barcodes(seq: str, 
                   ordered: bool = False) -> Generator[str]:

    """Generate a stream of barcodes encoding an amino
    acid sequence.

    Makes no consideration of codon usage preferences. If `ordered` is `True`,
    it is ignored if the number of possible combinations is more than 
    100,000.

    Parameters
    ----------
    seq : str
        Amino acid sequence to encode, in one-letter code.
    ordered : bool
        Whether to produce barcodes in sorted order. Default: False.

    Yields
    ------
    sequence : str
        DNA sequence encoding amino acid sequence.

    Examples
    --------
    >>> list(codon_barcodes("L", ordered=True))  # doctest: +NORMALIZE_WHITESPACE
    ['CTT', 'CTC', 'CTA', 'CTG', 'TTA', 'TTG']
    >>> list(codon_barcodes("L"))  # doctest: +SKIP
    ['TTA', 'CTT', 'CTA', 'CTG', 'CTC', 'TTG']

    """

    codons = [_CODONS[aa] for aa in seq]
    n_combos = reduce(operator.mul, map(len, codons))
    combos_tried = set()

    if ordered and n_combos < 1e5:
        
        codons = (_CODONS[aa] for aa in seq)
            
        for combo in product(*codons):

            yield ''.join(combo)

    else:

        while len(combos_tried) < n_combos:

            this_sample = ''.join(sample(codon, k=1)[0] for codon in codons)

            if this_sample not in combos_tried:

                combos_tried.add(this_sample)
                
                yield this_sample          


def infinite_barcodes(length: int = 12,
                      alphabet: Union[Iterable[str], Iterable[Mapping]] = sq.sequences.DNA,
                      check_used: bool = True) -> Generator[str]:

    """Generate an stream of random barcodes by 
    randomly sampling from an alphabet.

    Not actually infinite by default. Set `check_used = False`. This
    will produce barcodes forever, so make sure you have some 
    end condition in your loop.

    Parameters
    ----------
    length : int
        Length of barcode to generate.
    alphabet : Iterable, optional
        Set of letters from which to sample.
    check_used : bool
        Only produce unique sequences. Default: True.

    Yields
    ------
    sequence : str
        Sequence with desired length.

    Examples
    --------
    >>> sorted(infinite_barcodes(2))  # doctest: +NORMALIZE_WHITESPACE
    ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
    >>> sorted(infinite_barcodes(2, alphabet='cats'))  # doctest: +NORMALIZE_WHITESPACE
    ['aa', 'ac', 'as', 'at', 'ca', 'cc', 'cs', 'ct', 'sa', 'sc', 'ss', 'st', 'ta', 'tc', 'ts', 'tt']
    >>> sorted(infinite_barcodes(2, alphabet=transition_matrix(['ATCG', 'ATTT'])))  # doctest: +NORMALIZE_WHITESPACE
    ['AT']
    >>> for bc in infinite_barcodes(20, check_used=False):  # doctest: +SKIP
    ...     print(bc)
    ...     break
    ... 
    ATCAGTCGTCACACTAGTTA

    """
    
    combos_tried = set()

    if isinstance(alphabet[0], str):
        
        ones = (1., ) * len(alphabet)
        alphabet = (({None: (alphabet, ones)},) + 
                    tuple({letter: (alphabet, ones) for letter in alphabet} 
                           for _ in range(length - 1)))
        
    alphabet = alphabet[:length]
    n_combos = reduce(operator.mul, 
                      (max(len(letters) for _, (letters, _) in position.items()) 
                       for position in alphabet))

    while len(combos_tried) < n_combos or not check_used:

        this_sample = []
        letter = None

        for transition_matrix in alphabet:
            
            letters, weights = transition_matrix[letter]
                
            letter = choices(letters, weights=weights, k=1).pop()
            this_sample.append(letter)

        this_sample = ''.join(this_sample)

        if check_used and (this_sample not in combos_tried):

            combos_tried.add(this_sample)

            yield this_sample

        elif not check_used:

            yield this_sample
