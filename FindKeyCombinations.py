# -*- coding:utf-8 -*-
import PKIpipeline as PKI
import pandas as pd
import numpy as np
import random
# read file

networkstructure = pd.read_csv('PKI_submit\\ESCdataset\\Network.csv', header=0)
derivdata = pd.read_csv("PKI_submit\\ESCdataset\\ChangeRate.csv", index_col=0)  # col time
timeseriesdata = pd.read_csv("PKI_submit\\ESCdataset\\TimeSeriesExpressionData.csv", index_col=0)

print(networkstructure.head(), "\n", derivdata.head(), "\n", timeseriesdata.head())

# knockout list should be 0,1; length equal to TFlist
#### input: knockout_index, TFlist, relationshipuse, data, derivdata, sumscore_origin
#### output: delete,ratio_sum
relationship_input, TF = PKI.preprocess2(networkstructure)
sum_score0, mean_socre0 = PKI.Rsquared_origin(relationship=relationship_input,
                                                                                 data=timeseriesdata.T,
                                                                                 derivdata=derivdata.T)

pgm = 0.2  # possibility of gene mutation

pgc = 0.8  # possibility of crossover

breeding_algebra = 100

# Heuristic Initialization: add high ranked TF
a = [1] * 30
a[2] = 0
a[21] = 1
pop = [a for i in range(60)]


def objvaluefunc(pop, set_n):
    objvalue = []

    for i in range(len(pop)):

        x = pop[i]
        p = PKI.pki(x, TF, relationshipuse=relationship_input, data=timeseriesdata.T, derivdata=derivdata.T,
                sumscore_origin=sum_score0)[1]

        if np.sum(x) < 30 - set_n:
            p = p - 0.9

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

    return [bestvalue, bestindividual]


def sumfitvalue(fitvalue):
    total = 0

    for i in range(len(fitvalue)):
        total += fitvalue[i]

    return total


def cumsumfitvalue(fitvalue):
    t = 0

    for i in range(len(fitvalue)):
        t += fitvalue[i]
        fitvalue[i] = t


def selection(pop, fitvalue):  # 自然选择，轮盘赌算法

    totalvalue = sumfitvalue(fitvalue)

    newfitvalue = []

    for i in range(len(fitvalue)):
        m = fitvalue[i] / totalvalue
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

        if (ms[msin] < newfitvalue[fitin]):

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
    n = 2  # The number of TFs in combination. n=2,3,4
    all_result = []
    print("start")
    for i in range(0, 10):

        results = []
        for i in range(breeding_algebra):
            objvalue = objvaluefunc(pop, set_n=n)

            fitvalue = fitvaluefunc(objvalue)

            totalvalue = sumfitvalue(fitvalue)

            [bestvalue, bestindividual] = best(pop, fitvalue)

            results.append([bestvalue, bestindividual])

            # print([bestvalue, bestindividual])

            selection(pop, fitvalue)

            cross(pop, pgc)

            muta(pop, pgm)

        results.sort()

        print(results[-1])
        all_result.append(results[-1])

print("Number of TFs:", n, "\nKey Combination Result:", all_result)
# import pandas as pd
# a=pd.DataFrame(all_result)
# a.to_csv("results3.csv")
