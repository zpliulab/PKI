#-*- coding: utf-8 -*-

import pandas as pd


data = pd.read_csv("log1pdata.csv",index_col=0,header=0)

t0 = data.iloc[:, 0:91].mean(1)
t12 = data.iloc[:, 92:193].mean(1)
t24 = data.iloc[:, 194:259].mean(1)
t36 = data.iloc[:, 260:431].mean(1)
t72 = data.iloc[:, 432:569].mean(1)
t96 = data.iloc[:, 570:757].mean(1)

tf = pd.read_csv('setTF.csv', index_col=0,header=0)

timedata=pd.concat([t0,t12,t24,t36,t72,t96], axis=1)

sixpointdata=pd.DataFrame()

for i in tf['0']:
    if i in data._stat_axis.values.tolist():
        sixpointdata[i]=timedata.loc[i]
sixpointdata=sixpointdata.T
#sixpointdata.to_csv("sixpointdata.csv")

def diffdata(dataframe):
    data = pd.DataFrame()
    n=len(dataframe.ix[0])
    for i in range(0,n-1):
        data[i]=dataframe.iloc[:,i+1]-dataframe.iloc[:,i]
    data.to_csv("diffdata.csv")
    return data

#diffdata(sixpointdata)

def process_rel(filename): #rel_f
    relationship = pd.read_csv(filename)
    all_target=[]

    relationship_f=[]

    for i in list(relationship['target']):
        if i not in all_target:
            all_target.append(i)

    for i in all_target:
        one=[]
        one.append(i)
        for j in range(0, len(relationship['target'])):
            if i==relationship["target"][j]:
                one.append(relationship['tf'][j])
                print one
        relationship_f.append(one)

    print relationship_f
    relationship_f=pd.DataFrame(relationship_f)
    relationship_f.to_csv("relationship_f")

process_rel("rel_f1.csv")




