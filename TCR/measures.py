import numpy as np
import itertools as it
import operator as op
from distances import *

def error_evaluation(weights, data, gaps=0):
    N1, N2 = map(len, data)
    data = data.values()
    M1 = list(it.accumulate(
            blosum62_distance(data[0], data[0], weights, gaps),
            op.add))[-1]
    M2 = list(it.accumulate(
            blosum62_distance(data[0], data[1], weights, gaps),
            op.add))[-1]
    M3 = list(it.accumulate(
            blosum62_distance(data[1], data[1], weights, gaps),
            op.add))[-1]
    return N1*N2*(M1/N1**2+M3/N2**2)/M2
