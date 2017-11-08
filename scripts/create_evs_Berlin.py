#!/usr/bin/env python
import pandas as pd
import sys
import numpy as np

FILENAME=sys.argv[1]
OUTPUTNAME1=sys.argv[2]
OUTPUTNAME2=sys.argv[3]

data = pd.read_csv(FILENAME, sep = '\t', header=None, index_col = 0).T

EV1 = pd.DataFrame({'onset': data['left_1'], 'duration': 0, 'value': 1}) 
EV2 = pd.DataFrame({'onset': data['right_1'], 'duration': 0, 'value': 1}) 

EV1.to_csv(OUTPUTNAME1, header=False, columns=['onset', 'duration', 'value'], index=False, sep = ' ')
EV2.to_csv(OUTPUTNAME2, header=False, columns=['onset', 'duration', 'value'], index=False, sep = ' ')

