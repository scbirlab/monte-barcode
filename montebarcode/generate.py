"""Functions for generating random barcodes."""

from collections.abc import Generator, Iterable
from functools import reduce
from itertools import product
import operator
from random import choices, sample

import streq as sq

from .utils import _CODONS


def codon_barcodes(seq: str, ordered: bool = False) -> Generator[str]:

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
                      alphabet: Iterable[str] = sq.sequences.DNA,
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
    >>> sorted(infinite_barcodes(2))  # doctest: +SKIP
    ['AA', 'AG', 'AG', 'AT', 'CA', 'CA', 'CA', 'CC', 'CG', 'CT', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG']
    >>> for bc in infinite_barcodes(20, check_used=False):  # doctest: +SKIP
    ...     print(bc)
    ...     break
    ... 
    ATCAGTCGTCACACTAGTTA

    """

    n_combos = len(alphabet) ** length
    combos_tried = set()

    while len(combos_tried) < n_combos or not check_used:

        this_sample = ''.join(choices(alphabet, k=length))

        if not check_used or (this_sample not in combos_tried):

            combos_tried.add(this_sample)

            yield ''.join(choices(alphabet, k=length))
