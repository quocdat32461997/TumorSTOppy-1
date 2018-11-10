from importlib import reload
import numpy as np
import itertools as it
from tumorstoppy import distances, measures, data

cdr3 = data.CDR3_13
error_evaluation = measures.error_evaluation
#print(error_evaluation(None, cdr3['training'], verbose=True))
