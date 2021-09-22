# -*- coding: utf-8 -*
# ------------------------------------------------------------------------------
# python 2 and 3 compatible
from __future__ import division
import numpy as np
# ------------------------------------------------------------------------------



def autocorr(x):
    result = np.correlate(x, x, mode='full')
    return result[result.size // 2:]
