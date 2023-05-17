import doctest
import montebarcode as mb

if __name__ == '__main__':

    doctest.testmod(mb.checks)
    doctest.testmod(mb.generate)