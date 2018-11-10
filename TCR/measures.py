import numpy as np
import itertools as it
import operator as op
from distances import *

def error_evaluation(weights, data, gaps=0):
    N1, N2 = map(len, data)
    M1, M2, M3 = 0, 0, 0
    data = data.values()

    for i in range(N1):
        for j in range(i+1,N1):
            M1 += blosum62_distance(tumor[i],tumor[j],weights)
    M1 /= N1*N1

    for i in range(N1):
        for j in range(N2):
            M2 += blosum62_distance(tumor[i],nontumor[j],weights)
    M2 /= N1*N2

    for i in range(N2):
        for j in range(i+1,N2):
            M3 += distance(nontumor[i],nontumor[j],weights)
    M3 /= N2*N2

    return (M1 + M3) / M2
