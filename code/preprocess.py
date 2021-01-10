# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from biolearns import preprocessing

workdir = '' # change to your directory

mRNA = pd.read_excel(workdir+'73CN-AML-RNA-TCGA.xls', index_col=0)
miRNA = pd.read_excel(workdir+'73CN-AML-miRNA-TCGA.xls', index_col=0)
meth = pd.read_csv(workdir+'73CN-AML-cn_meth-TCGA.csv', index_col=0)

mRNA.columns = [int(c.split('-')[2]) for c in mRNA.columns]
miRNA.columns = [int(c.split('-')[2]) for c in miRNA.columns]
meth.columns = [int(c) for c in meth.columns]

mRNA.index = ['mRNA*' + i for i in mRNA.index]
miRNA.index = ['miRNA*' + i for i in miRNA.index]
meth.index = ['meth*' + i for i in meth.index]

mRNA = preprocessing.expression_filter(mRNA, 0.2, 0.2)
miRNA = preprocessing.expression_filter(miRNA, 0.2, 0.2)
meth = preprocessing.expression_filter(meth, 0.8, 0.8)

joint = mRNA.append(miRNA) 
joint = joint.append(meth) 

joint.to_csv(workdir+'joint.csv')

joint_normalized = (joint.T / np.linalg.norm(joint, axis=1)).T
joint_normalized.to_csv(workdir+'joint_normalized.csv')
