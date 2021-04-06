# -*- coding:utf-8 -*-
import pandas as pd
import numpy as np


dataall=pd.read_csv("GSE75748_sc_time_course_ec.csv",index_col=0,header=0).apply(np.log1p)
dataall.to_csv("log1pdata.csv")
tf = pd.read_table('tf_new.txt', index_col=None,header=None)
data=pd.DataFrame()
for i in tf[0]:
    if i in dataall._stat_axis.values.tolist():
        data[i]=dataall.loc[i]
data=data.T
# data=pd.read_csv("D:\\A\work\\2019\\scrna\\data3\\data.txt",header=None,sep="\t").T
# tf = pd.read_table('tf.txt', index_col=None,header=None)
#
# data.columns=tf[0]
# data=data.T
#data=data.replace(0,np.nan)
t0 = data.iloc[:, 0:91].mean(1)
t12 = data.iloc[:, 92:193].mean(1)
t24 = data.iloc[:, 194:259].mean(1)
t36 = data.iloc[:, 260:431].mean(1)
t72 = data.iloc[:, 432:569].mean(1)
t96 = data.iloc[:, 570:757].mean(1)
timedata=pd.concat([t0,t12,t24,t36,t72,t96], axis=1)
timedata.to_csv("sixpoint.csv")
timedata['variances']=timedata.var(axis = 1)
scodedataall=timedata.sort_values("variances",ascending=False)
tf200= pd.DataFrame(scodedataall.iloc[0:200,6])
tf100=pd.read_csv("tf.txt",header=None)
tf200name=tf200._stat_axis.values.tolist()
print (tf200name)
p=0
for i in tf100[0]:
    if i in tf200name:
        p=p+1
        print (i)
        print (p)
tf200.to_csv("200TF.csv")

