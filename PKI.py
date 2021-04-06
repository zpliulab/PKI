# -*- coding:utf-8 -*-

import random
import numpy as np
from sklearn.linear_model import RidgeCV
import pandas as pd


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

print (sum_score0)

pgm = 0.2  # possibility of gene mutation

pgc = 0.8  # possibility of crossover

breeding_algebra = 100

#Heuristic Initialization
a=[1]*97
a[60]=0
a[93]=0
pop = [a for i in range(60)]


#fitness fuction
def rsquared(V):

    delete = []
    delete_row = []
    for p in range(0, len(tfs['0'])):

        if V[p] == 0:

            delete.append(tfs.iloc[p,0])
            delete_row.append(tfs.iloc[p,0])

    diffrelationship = pd.read_csv("relationship_f", header=0, index_col=0)

    mark=[]
    for i in range(0,len(diffrelationship)):

        if(len(delete_row)>0):

            for w in range(0,len(delete_row)):

                if delete_row[w]==diffrelationship.iloc[i,0]:

                    mark.append(i)
                    break
    diffrelationship = diffrelationship.drop(index=mark)

    for i in range(0,len(diffrelationship)):

        for j in range(0,len(diffrelationship.iloc[i,:])):

            if(len(delete)>0):

                for k in range(0,len(delete)):

                    if delete[k]==diffrelationship.iloc[i,j]:
                        diffrelationship.iloc[i,j] = np.NaN


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

    return (sum_score0 - sum_score)/sum_score0


def objvaluefunc(pop,set_n):
    objvalue = []

    for i in range(len(pop)):

        x = pop[i]
        p=rsquared(x)

        if np.sum(x)<97-set_n:
            p=p-0.9

        objvalue.append(p)

    return objvalue


def fitvaluefunc(objvalue):


    fitvalue = []

    for i in range(len(objvalue)):

        if (objvalue[i] > 0):

            fitvalue.append(objvalue[i])

        else:

            fitvalue.append(0)

    return fitvalue

def best(pop, fitvalue):  # 计算最适合个体和最适应值

    bestvalue = fitvalue[0]
    bestindividual = pop[0]

    for i in range(1, len(fitvalue)):

        if (fitvalue[i] > bestvalue):

            bestvalue = fitvalue[i]
            bestindividual = pop[i]

    return [bestvalue,bestindividual]

def sumfitvalue(fitvalue):

    total = 0

    for i in range(len(fitvalue)):

        total += fitvalue[i]

    return total

def cumsumfitvalue(fitvalue):

    t=0

    for i in range(len(fitvalue)):

        t+=fitvalue[i]
        fitvalue[i]=t

def selection(pop, fitvalue):  # 自然选择，轮盘赌算法

    totalvalue = sumfitvalue(fitvalue)

    newfitvalue = []

    for i in range(len(fitvalue)):

        m=fitvalue[i]/totalvalue
        newfitvalue.append(m)

    cumsumfitvalue(newfitvalue)

    ms = []

    poplen = len(pop)

    for i in range(poplen):

        ms.append(random.random())

    ms.sort()

    msin = 0
    fitin = 0
    newpop = pop

    while msin < poplen:

        if(ms[msin] < newfitvalue[fitin]):

            newpop[msin] = pop[fitin]
            msin = msin + 1

        else:

            fitin = fitin + 1

    pop = newpop

def cross(pop, pgc):  # crossover

    for i in range(len(pop) - 1):

        if (random.random() < pgc):

            cnum = random.randint(0, len(pop[i]))

            t1 = []

            t2 = []

            t1.extend(pop[i][0:cnum])

            t1.extend(pop[i + 1][cnum:len(pop[i])])

            t2.extend(pop[i + 1][0:cnum])

            t2.extend(pop[i][cnum:len(pop[i])])

            pop[i] = t1

            pop[i + 1] = t2

def muta(pop, pgm):  # gene mutation

    for i in range(len(pop)):

        if (random.random() < pgm):

            mnum = random.randint(0, len(pop[i]) - 1)

            if (pop[i][mnum] == 0):

                pop[i][mnum] = 1

            else:

                pop[i][mnum] = 0


if __name__ == "__main__":
    n=4  # The number of TFs in combination. n=1,2,3,4
    all_result=[]

    for i in range(0,100):
        results = []
        for i in range(breeding_algebra):
            objvalue = objvaluefunc(pop,set_n=n)

            fitvalue = fitvaluefunc(objvalue)

            totalvalue = sumfitvalue(fitvalue)

            [bestvalue, bestindividual] = best(pop, fitvalue)

            results.append([bestvalue, bestindividual])

            print([bestvalue, bestindividual])

            selection(pop, fitvalue)

            cross(pop, pgc)

            muta(pop, pgm)

        results.sort()

        print(results[-1])
        all_result.append(results[-1])


print(all_result)
import pandas as pd
a=pd.DataFrame(all_result)
a.to_csv("results4.csv")


