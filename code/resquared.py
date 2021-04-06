# -*- coding:utf-8 -*-


from sklearn.linear_model import RidgeCV
import pandas as pd
import numpy as np

def _read_file(filename):
    lists=[]
    f = open(filename, "r")
    list_raw = f.readlines()
    for fields in list_raw:
        fields = fields.strip()
        fields = fields.split(" ")
        lists.append(fields)
    return lists

def combinations(iterable, r):
    # combinations('ABCD', 2) --> AB AC AD BC BD CD
    # combinations(range(4), 3) --> 012 013 023 123
    pool = tuple(iterable)
    n = len(pool)
    if r > n:
        return
    indices = list(range(r))
    yield tuple(pool[i] for i in indices)
    while True:
        for i in reversed(range(r)):
            if indices[i] != i + n - r:
                break
        else:
            return
        indices[i] += 1
        for j in range(i+1, r):
            indices[j] = indices[j-1] + 1
        yield tuple(pool[i] for i in indices)


diffrelationship=pd.read_csv("relationship_f",header=0,index_col=0)
tfs=pd.read_csv("tf.csv",header=0,index_col=0)
targetsname=diffrelationship['0']
diffdata=pd.read_csv("diffdata.csv",header=0,index_col=0)
mydata=pd.read_csv("sixpointdata.csv",header=0,index_col=0)

mydata = mydata.T
mydata = mydata.drop(['5'])
diffdata = diffdata.T

constant_diffrelationship=pd.read_csv("relationship_f",header=0,index_col=0)

score_list = []

for i in range(0, len(diffrelationship)):
    y_name = diffrelationship.iloc[i,0]
    x_name = diffrelationship.iloc[i,1:].dropna()
    y = diffdata[y_name]
    X = mydata[x_name]
    clf = RidgeCV(fit_intercept=False)
    clf.fit(X, y)
    score = clf.score(X, y)
    score_list.append(score)
sum_score0 = np.mean(score_list)


def rsquared(V):
    print("111111111111111111111111")

    delete = []
    delete_row = []
    for p in range(0, len(tfs['0'])):
        if V[p] == 0:
            print('yes!')
            delete.append(tfs.iloc[p,0])
            print(delete)
            delete_row.append(tfs.iloc[p,0])

    diffrelationship = pd.read_csv("relationship_f", header=0, index_col=0)

    mark=[]
    for i in range(0,len(diffrelationship)):
        if(len(delete_row)>0):
            for w in range(0,len(delete_row)):
                if delete_row[w]==diffrelationship.iloc[i,0]:
                    #print("delete row")
                    mark.append(i)
                    break
    diffrelationship = diffrelationship.drop(index=mark)
    for i in range(0,len(diffrelationship)):

        for j in range(0,len(diffrelationship.iloc[i,:])):
            # print("start",diffrelationship[i])

            if(len(delete)>0):
                for k in range(0,len(delete)):
                    # print("k",k)
                    #print(delete[k])
                    #print(diffrelationship.iloc[i,j])
                    if delete[k]==diffrelationship.iloc[i,j]:
                        diffrelationship.iloc[i,j] = np.NaN
                        #print("delete")


    score_list = []

    for q in range(0, len(diffrelationship)):
        if len(diffrelationship.iloc[q,:].dropna()) < 2:
            continue
        else:
            y_name = diffrelationship.iloc[q,0]
            x_name = diffrelationship.iloc[q,1:].dropna()
            y = diffdata[y_name]
            X = mydata[x_name]
            clf = RidgeCV(fit_intercept=False)
            clf.fit(X, y)
            score = clf.score(X, y)
            score_list.append(score)

    sum_score = np.mean(score_list)
    num= (sum_score0 - sum_score)/sum_score0
    return num




#################################################
# resultlist=[]
# for i in range(0,len(tfs['0'])):
#     list = [1] * len(tfs['0'])
#     list[i]=0
#     print rsquared(list)
#     resultlist.append(rsquared(list))
# resultlist=pd.DataFrame(resultlist)
# a=pd.DataFrame()
# a['name']=tfs['0']
# a['num']=resultlist
# a=a.sort_values("num")
# a.to_csv("tfsort")

#################################################
# two_resultlist=[]
# outlist=[]
# for i in range(0,len(tfs["0"])-1):
#     for j in range(i+1,len(tfs["0"])):
#         outlist.append([i,j])
#     print outlist
#
#
# for i in outlist:
#     tflist = [1] * len(tfs["0"])
#     x1=i[0]
#     x2=i[1]
#     tflist[x1]=0
#     tflist[x2]=0
#     two_resultlist.append(rsquared(tflist))
# two_resultlist=pd.DataFrame(two_resultlist)
# two_resultlist.to_csv("tf2result")
# two_resultlist["comb2"]=outlist
# two_resultlist=two_resultlist.sort_values['0']
# two_resultlist.to_csv("sortTF2")



#three_resultlist=[]
# for i in list(combinations(range(len(tfs["0"])),3)):
#     tflist = [1] * len(tfs["0"])
#     x1=i[0]
#     x2=i[1]
#     x3=i[2]
#     tflist[x1]=0
#     tflist[x2]=0
#     tflist[x3] = 0
#     three_resultlist.append(rsquared(tflist))
# a=pd.DataFrame()
# a["num"]=three_resultlist
# a.to_csv("tf3result")
# a["comb3"]=list(combinations(range(len(tfs["0"])),3))
# a=a.sort_values(['num'],ascending=False)
# a.to_csv("sortTF3")
















