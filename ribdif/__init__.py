# -*- coding: utf-8 -*-
"""RibDif2 - Ribosomal differentiation v2

Documentation: https://github.com/Rob-murphys/RibDif2

"""
__authors__ = "Robert Murphy"
__licence__ = "GPL-3.0"
__version__ = (2, 0, 0)


# Try on python v3.5
import sys as _sys
if _sys.version_info[:2] < (3, 7):
    raise ImportError('Python version must be >= 3.7')