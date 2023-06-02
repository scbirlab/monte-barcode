"""Functions to check that barcodes conform."""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping, Sequence
from collections import Counter, defaultdict
from itertools import zip_longest
import sys

import streq as sq

from .utils import pprint_dict


def _print_rejection_reasons(x: Mapping, 
                            n_tried: str) -> None:
    
    key_val_str = {name: (s / n_tried) 
                   for name, s in x.items()}

    pprint_dict(key_val_str, '\nRejection reasons')

    return None


def minmax_distance(x: Sequence[str],
                    use_levenshtein: bool = True) -> Sequence[int]:
    
    """Get minimum and maximum distances among a set of barcodes.

    This compares all the pairwise distances except for self-self (diagonal
    of the distance matrix). If there are repeated barcodes in the list, 
    then the minimum distance will be zero.

    Parameters
    ----------
    x : list | tuple | set
        Set of barcodes to check.
    use_levenshtein : bool, optional
        Whether to use the Levenshtein distance. Default: True.

    Returns
    -------
    tuple
        Two integers: the first is the minimum distance, the second is the maximum distance.

    Examples
    --------
    >>> minmax_distance(['AAA', 'AAA'])
    (0, 0)
    >>> minmax_distance(['AAA', 'TCG', 'AAT'])
    (1, 3)
    >>> minmax_distance(['AAA', 'TCG', 'AAAT'], use_levenshtein=False)
    (0, 3)
    >>> minmax_distance(['AAA', 'TCG', 'AAAT'])
    (1, 4)

    """
    
    dist_function = sq.levenshtein if use_levenshtein else sq.hamming

    all_distances = [dist_function(seq1, seq2) 
                     for i, seq1 in enumerate(x) 
                     for j, seq2 in enumerate(x) if i != j]
    
    try:
        return min(all_distances), max(all_distances)
    except ValueError:
        return 0, 0
        

def Distance(min_distance: int = 2,
             use_levenshtein: bool = True) -> Callable[[str, Iterable], bool]:
    
    """Create a distance checking function.

    Uses the parameters to produce a function which takes as
    arguments a candidate barcode and a working list of barcodes,
    and returns True if the parameter conditions are not met.

    Parameters
    ----------
    min_distance : int
        The minimum distance allowed among all barcodes.
    use_levenshtein : bool, optional
        Whether to use the Levenshtein distance. Default: True.

    Returns
    -------
    function
        Checking function.

    Examples
    --------
    >>> Distance(1)('ATA', ['TCG', 'AAT'])
    False
    >>> Distance(2)('AAA', ['TCG', 'AAT'])
    True

    """
    
    dist_function = sq.levenshtein if use_levenshtein else sq.hamming

    def distance(word: str, 
                 corpus: Iterable) -> bool:

        for ref in corpus:
            
            hamming_dist = dist_function(word, ref)
            
            if hamming_dist < min_distance:

                return True
        
        return False

    return distance


def Identities() -> Callable[[str, Iterable], bool]:

    """Create an identity checking function.

    Produces a function which takes as arguments a candidate barcode 
    and a working list of barcodes, and returns True if the candidate 
    or its reverse complement is in the working list.

    Returns
    -------
    function
        Checking function.

    Examples
    --------
    >>> Identities()('AAA', ['TCG', 'AAT'])
    False
    >>> Identities()('AAA', ['TCG', 'AAA'])
    True
    >>> Identities()('AAA', ['TCG', 'TTT'])
    True

    """

    def identity(word: str, 
                 corpus: Iterable) -> bool:

        return word in corpus or sq.reverse_complement(word) in corpus
    
    return identity


def Palindrome() -> Callable[[str, Iterable], bool]:

    """Create a palindrome checking function.

    Produces a function which takes as arguments a candidate barcode 
    and a working list of barcodes, and returns True if the candidate 
    is palindromic. (Working list is ignored.)

    Returns
    -------
    function
        Checking function.

    Examples
    --------
    >>> Palindrome()('AAA', [])
    False
    >>> Palindrome()('AATT', [])
    True

    """

    def palindrome(word: str, _) -> bool:

        return word == sq.reverse_complement(word)

    return palindrome


def GCcontent(min: float = .35, 
              max: float = .65) -> Callable[[str, Iterable], bool]:
    
    """Create GC content checking function.

    Produces a function which takes as arguments a candidate barcode 
    and a working list of barcodes, and returns True if the candidate 
    is not within the bounds. (Working list is ignored.)

    Parameters
    ----------
    min : float
        Minimum acceptable proportion of GC content.
    max : float
        Maximum acceptable proportion of GC content.

    Returns
    -------
    function
        Checking function.

    Examples
    --------
    >>> GCcontent()('AATT', [])
    True
    >>> GCcontent()('AACG', [])
    False
    >>> GCcontent()('GGCG', [])
    True

    """

    def gc_content(word: str, _) -> bool:

        gc_ = sum([letter in 'GC' for letter in word]) / len(word)

        return gc_ < min or gc_ > max

    return gc_content


def _illumina_color_balance(x: Iterable[str], **kwargs) -> bool:

    (image1_avg, 
     image2_avg, 
     dark_avg) = ([sum(b in channel for b in bases) / len(bases) 
                   for bases in zip(*x)] for channel in dict(**kwargs).values())
    
    both_images = any(a <= 0. or b <= 0. for a, b in zip(image1_avg, image2_avg))
    signal = any(a >= 1. for a in dark_avg)
    
    return both_images or signal


def base_usage(x: Iterable[str]) -> Mapping:

    """Calculate the proportional base usage for a set of barcodes.
    
    Returns a dictionary mapping position along sequence to
    the distribution among the bases.

    Does not check for legitimate DNA alphabet.

    Parameters
    ----------
    x : Iterable
        Set of barcodes to check.

    Returns
    -------
    dict
        Base usage per position.

    Examples
    --------
    >>> base_usage(['AAA', 'TTT', 'GCT', 'CCA'])[0]['A']
    0.25
    >>> base_usage(['AAA', 'TTT', 'GCT', 'CCA'])[1]['G']
    0
    >>> base_usage(['AAA', 'TTT', 'GCT', 'CCA'])[2]['A']
    0.5

    """

    counter = defaultdict(Counter)

    for i, bases in enumerate(zip_longest(*x)):
        
        for j, base in enumerate(bases):

            counter[i][base] += 1

        for base in counter[i]:

            counter[i][base] /= (j + 1)

    return counter


def IlluminaColorBalance(green_4ch='GT', 
                         red_4ch='AC', 
                         green_2ch='AT', 
                         red_blue_2ch='AC',
                         image1_1ch='AT', 
                         image2_1ch='C',
                         alphabet=sq.sequences.DNA) -> Callable[[str, Iterable], bool]:
    
    """Create a color balance checking function.

    Produces a function which takes as arguments a candidate 
    barcode and a working list of barcodes, and returns True 
    if adding the barcode to the working set would give suboptimal
    color balance across a range of Illumina SBS platforms.

    This is a very stringent check in order to ensure barcodes
    will wrok on multiple platforms.

    If there are 4 or fewer barcodes in total, the checking function 
    will only check for two dark bases (for 1 channel chemistry) at the 5' end.

    If there are at least 5 barcodes in total, the checking function
    will test for color balance among the channels and images.

    Parameters
    ----------
    green_4ch : str, optional
        Bases detected by the green channel in 4 channel chemistry. Default: 'GT', 
    red_4ch : str, optional
        Bases detected by the red channel in 4 channel chemistry. Default: 'AC', 
    green_2ch : str, optional
        Bases detected by the green channel in 2 channel chemistry. Default: 'AT', 
    red_blue_2ch : str, optional
        Bases detected by the red/blue channel in 2 channel chemistry. Default: 'AC',
    image1_1ch : str, optional
        Bases detected by the first image in 1 channel chemistry. Default: 'AT', 
    image2_1ch : str, optional
        Bases detected by the second image in 1 channel chemistry. Default: 'C',
    alphabet : str, optional
        Letters comprising the alphabet. Used to infer the dark bases for each chemistry.

    Returns
    -------
    function
        Checking function.

    Examples
    --------
    >>> IlluminaColorBalance()('GGAT', ['TCGC', 'AAAG'])
    True
    >>> IlluminaColorBalance()('AAAT', ['TCGC', 'AAAG'])
    False
    >>> IlluminaColorBalance()('AAAT', ['TCGC', 'ACAG', 'TGGC', 'ATCG'])
    True
    >>> IlluminaColorBalance()('AAAT', ['TCGC', 'CCAG', 'TGGC', 'ATCG'])
    False

    """

    dark_4ch  = ''.join(letter for letter in alphabet 
                        if letter not in (green_4ch + red_4ch))
    dark_2ch = ''.join(letter for letter in alphabet 
                        if letter not in (green_2ch + red_blue_2ch))
    dark_1ch  = ''.join(letter for letter in alphabet 
                        if letter not in (image1_1ch + image2_1ch))
    dark_x2 = dark_1ch * 2

    def color_balance(word: str, corpus: Iterable) -> bool:

        start_GG = word.startswith(dark_x2)

        proposed_set = corpus + [word]

        if len(proposed_set) < 5:

            return start_GG
        
        else:

            _4ch = _illumina_color_balance(proposed_set, 
                                           green=green_4ch, 
                                           red=red_4ch, 
                                           dark=dark_4ch)
            _2ch = _illumina_color_balance(proposed_set, 
                                           green=green_2ch, 
                                           red_blue=red_blue_2ch, 
                                           dark=dark_2ch)
            _1ch = _illumina_color_balance(proposed_set, 
                                           image1=image1_1ch, 
                                           image2=image2_1ch, 
                                           dark=dark_1ch)

            return any((start_GG, _4ch, _2ch, _1ch))

    return color_balance


def Homopolymer(length: int = 4) -> Callable[[str, Iterable], bool]:

    """Create a homopolymer checking function.

    Produces a function which takes as arguments a candidate barcode 
    and a working list of barcodes, and returns True if the candidate 
    contains a homopolymer `length` or longer. (Working list is ignored.)

    Parameters
    ----------
    length : int
        Minimum length of homopolymer to detect.

    Returns
    -------
    function
        Checking function.

    Examples
    --------
    >>> Homopolymer(3)('AAAT', [])
    True
    >>> Homopolymer(4)('AAAT', [])
    False

    """

    hps = set([letter * length for letter in 
               sq.sequences.DNA])

    def homopolymer(word: str, _) -> bool:

        return any(hp in word for hp in hps)

    return homopolymer


def RestrictionSites(n: int = 0) -> Callable[[str, Iterable], bool]:

    """Create a Type IIS restriction site checking function.

    Produces a function which takes as arguments a candidate barcode 
    and a working list of barcodes, and returns True if the candidate 
    contains a Type IIS resitriction site commonly used in Golden
    Gate cloning. (Working list is ignored.)

    Parameters
    ----------
    n : int
        Minimum number of restriction sites to tolerate.

    Returns
    -------
    function
        Checking function.

    Examples
    --------
    >>> RestrictionSites()('AAATGGTCTC', [])
    True
    >>> RestrictionSites(1)('AAATGGTCTC', [])
    False
    >>> RestrictionSites()('AAATGCTCTC', [])
    False
    
    """
    
    def restriction_sites(word: str, _) -> bool:

        return sq.count_re_sites(word) > n

    return restriction_sites


def make_checks(barcodes: Iterable[str],
                n: int, 
                checks: Iterable[Callable[[str, Sequence], bool]],
                max_rejection_rate: float = 1.,
                max_tries: int = 10000,
                quiet: bool = False) -> Sequence[dict, int, Sequence[str]]:
    
    """Check barcode list conforms to the checks.

    Keeps a tally of rejection rate and rejection reasons.

    Parameters
    ----------
    barcodes : Iterable[str]
        List or generator of barcodes to check.
    n : int
        Minimum number of barcodes to accept. Stops when this is reached.
    checks : Iterable
        List of functions which take a candidate barcode and 
        previous barcodes as arguments and return True if the
        candidate should be rejected.
    length : int
        Length of barcodes to generate.
    max_rejection_rate : float, optional
        Rejection rate above which to terminate. Default: 1.
    max_tries : int, optional
        Number of barcodes to try before enforcing `max_rejection_rate`. Default 100000.
    quiet : bool, optional
        Whether to report progress. Default: True.
    
    Returns
    -------
    dict
        Counts of rejections based on failing each check.
    int
        Number of random barcode candidates tried.
    list
        List of n random barcodes passing the checks.

    Raises
    ------
    ValueError
        When rejection rate goes above the threshold.

    Examples
    --------
    >>> make_checks(['AAAT', 'TCGC', 'ACAG', 'TGGC', 'ATCG'], 5, checks=[IlluminaColorBalance()], quiet=True)  # doctest: +NORMALIZE_WHITESPACE
    (Counter({'color_balance': 1}), 5, ['AAAT', 'TCGC', 'ACAG', 'TGGC'])
    >>> make_checks(['AAAT', 'TCGC', 'CCAG', 'TGGC', 'ATCG'], 5, checks=[IlluminaColorBalance()], quiet=True)  # doctest: +NORMALIZE_WHITESPACE
    (Counter(), 5, ['AAAT', 'TCGC', 'CCAG', 'TGGC', 'ATCG'])
    >>> make_checks(['AAAT', 'TCGC', 'ACAG', 'TGGC', 'ATCG'], 4, checks=[IlluminaColorBalance()], quiet=True)  # doctest: +NORMALIZE_WHITESPACE
    (Counter(), 4, ['AAAT', 'TCGC', 'ACAG', 'TGGC'])
    >>> # 
    >>> checks = [Homopolymer(), Palindrome()]
    >>> make_checks(['AAAAT', 'CCCGGG', 'ATCGCG', 'GCCGAT'], 4, checks=checks, quiet=True)  # doctest: +NORMALIZE_WHITESPACE
    (Counter({'homopolymer': 1, 'palindrome': 1}), 4, ['ATCGCG', 'GCCGAT'])
    >>> make_checks(['AAAAT', 'CCCGGG', 'ATCGCG', 'GCCGAT'], 1, checks=checks, quiet=True)  # doctest: +NORMALIZE_WHITESPACE
    (Counter({'homopolymer': 1, 'palindrome': 1}), 3, ['ATCGCG'])

    """
    
    accepted_barcodes = []
    fail_tally = Counter()
    n_tried, n_rejected, n_accepted = 0, 0, 0
    rejection_rate = 0.

    for i, new_barcode in enumerate(barcodes):

        n_tried = i + 1

        check_fails = [check.__name__ for check in checks
                       if check(new_barcode, accepted_barcodes)]

        if any(check_fails):

            fail_tally += Counter(check_fails)
            n_rejected += 1

            if n_tried > max_tries and rejection_rate >= max_rejection_rate:
                
                if not quiet:
                    _print_rejection_reasons(fail_tally, n_tried)

                message = (f'Rejection rate reached {rejection_rate:.2f} '
                          f'({n_rejected} / {n_tried}), '
                          'which is above the maximum '
                          f'threshold of {max_rejection_rate:.2f}')
                
                raise ValueError(message)
            
        else:
            
            accepted_barcodes.append(new_barcode)
            n_accepted = len(accepted_barcodes)

        rejection_rate =  n_rejected / n_tried
        
        if not quiet:
            print(f'\r> Tried {n_tried} barcodes,',
                f'rejected {n_rejected}, accepted {n_accepted};',
                f'rejection rate is {rejection_rate:.2f}',
                file=sys.stderr, end='')
        
        if n_accepted >= n:

            if not quiet:
                print('', file=sys.stderr)
            break
    
    if not quiet:
        _print_rejection_reasons(fail_tally, n_tried)
        
    return fail_tally, n_tried, accepted_barcodes