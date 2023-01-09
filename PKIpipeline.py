# -*- coding: utf-8 -*-

"""
Created on 15/10/2022 PM10:09
@author: WANG Yijuan
@email: wangyijuan[at]mail.sdu.edu.cn
"""

# use dream to test model
# imput file
# 1.network structure (header=None, sep\t)
# 2.derive data (csv)
# 3.timeseries data (index_col=0, sep\t)

import pandas as pd
from sklearn.linear_model import RidgeCV
from sklearn.linear_model import LinearRegression as LR
import numpy as np


def preprocess(networkstructure, timeseriesdata, derivdata):
    # preprocess
    # 1.network structure to model relationship #tf target 1/0
    # 2.timeseriesdata replicant
    # 3.derivmatrix T, add columns
    # 4.generate TF list

    rel = networkstructure[networkstructure.iloc[:, 2] == 1]
    rel.columns = ["TF", "Target", "mark"]
    # print(rel)
    all_target = []
    relationship = []

    for i in list(rel['Target']):
        if i not in all_target:
            all_target.append(i)
    for i in all_target:
        one = []
        one.append(i)
        for j in range(0, len(rel['Target'])):
            if i == rel["Target"][j]:
                one.append(rel['TF'][j])
                # print(one)
        relationship.append(one)
    relationship = pd.DataFrame(relationship)
    # print(relationship)

    data2 = timeseriesdata.T
    timepoint = np.unique(data2.columns)
    data = pd.DataFrame()
    for i in timepoint:
        data[i] = data2[i].mean(1)
    data = data.T
    # print("data",data,data.columns)

    derivdata = derivdata.T
    derivdata.columns = data.columns

    TF = np.unique(rel["TF"])

    print("Preprocess done.")
    return relationship, data, derivdata, TF


def preprocess2(networkstructure):
    # preprocess
    # 1.network structure to model relationship #tf target 1/0
    # 2.timeseriesdata replicant
    # 3.derivmatrix T, add columns
    # 4.generate TF list
    rel = networkstructure
    rel.columns = ["TF", "Target"]
    # print(rel)
    all_target = []
    relationship = []

    for i in list(rel['Target']):
        if i not in all_target:
            all_target.append(i)
    for i in all_target:
        one = []
        one.append(i)
        for j in range(0, len(rel['Target'])):
            if i == rel["Target"][j]:
                one.append(rel['TF'][j])
                # print(one)
        relationship.append(one)
    relationship = pd.DataFrame(relationship)
    # print(relationship)
    TF = np.unique(rel["TF"])

    return relationship, TF


def Rsquared_origin(relationship, data, derivdata):
    # Before pseudo knockout
    score_list = []
    coeflist = []
    for i in range(0, len(relationship)):
        y_name = relationship.iloc[i, 0]
        x_name = relationship.iloc[i, 1:].dropna()
        y = derivdata[y_name]
        # print(y)
        X = data[x_name]
        # print(X)
        if len(relationship.iloc[i, :].dropna()) == 2:
            # print("choose lm")
            clf = LR(fit_intercept=True)
            clf.fit(X, y)
            score = clf.score(X, y)
            score_list.append(score)
            coeflist.append(clf.coef_)
            # print('score', score)
        else:
            # print("choose Ridge")
            clf = RidgeCV(fit_intercept=True, alphas=[0.00001, 0.005, 0.01, 0.05, 0.1])
            clf.fit(X, y)
            score = clf.score(X, y)
            score_list.append(score)
            coeflist.append(clf.coef_)
            # print('coef_matrix:\n', clf.coef_)
            # print('model:\n', clf)
            # print('alpha',clf.alpha_ )
            # print('score', score)
            # print('adscore', adscore)

    sum_score0 = np.sum(score_list)
    mean_score0 = np.mean(score_list)

    score_list = pd.DataFrame(score_list)
    score_list.to_csv("scorelist.csv")

    print("Origin R-squared calculation done.")
    # coeflist=pd.DataFrame(coeflist)
    # coeflist.to_csv("coeflist.csv")
    return sum_score0, mean_score0


# knockout list should be 0,1; length equal to TF list
#### input: knockout_index, TFlist, relationshipuse, data, derivdata, sumscore_origin
#### output: delete,ratio_sum
def pki(knockout_index, TFlist, relationshipuse, data, derivdata, sumscore_origin):
    print("enter pki")

    ###########pesudo knockout process#############
    relationshipuse_copy = relationshipuse.copy()
    delete = []
    delete_row = []
    for p in range(0, len(TFlist)):
        if knockout_index[p] == 0:
            print('delete No.', p)
            delete.append(TFlist[p])
            # print(delete)
            delete_row.append(TFlist[p])
    #########as target also delete#########
    mark = []
    for i in range(0, len(relationshipuse_copy)):
        if (len(delete_row) > 0):
            for w in range(0, len(delete_row)):
                if delete_row[w] == relationshipuse_copy.iloc[i, 0]:
                    # print("delete row",mark,len(mark))
                    mark.append(i)
                    break
    if (mark):
        print("enter delete row", mark)
        relationshipuse_copy = relationshipuse_copy.drop(index=mark)

    for i in range(0, len(relationshipuse_copy)):

        for j in range(1, len(relationshipuse_copy.iloc[i, :])):
            # when knockout gene as Target, RETAIN; or from 0
            # print("start",relationship[i])

            if (len(delete) > 0):
                for k in range(0, len(delete)):
                    # print("k",k)
                    # print(delete[k])
                    # print(relationship.iloc[i,j])
                    if delete[k] == relationshipuse_copy.iloc[i, j]:
                        relationshipuse_copy.iloc[i, j] = np.NaN
                        # print("delete")

    ########### modeling #############
    score_list = []
    # mergscore_list = []

    for q in range(0, len(relationshipuse_copy)):
        # print("here:", relationship.iloc[q,:])
        if len(relationshipuse_copy.iloc[q, :].dropna()) < 2:
            continue
        else:
            y_name = relationshipuse_copy.iloc[q, 0]
            x_name = relationshipuse_copy.iloc[q, 1:].dropna()
            # print(y_name)
            y = derivdata[y_name]
            X = data[x_name]
            if len(relationshipuse_copy.iloc[q, :].dropna()) == 2:
                # print("choose lm")
                clf = LR(fit_intercept=True)
                clf.fit(X, y)
                score = clf.score(X, y)
                score_list.append(score)
                # mergscore_list.append(score)

            else:
                # print("choose Ridge")
                clf = RidgeCV(fit_intercept=True, alphas=[0.00001, 0.005, 0.01, 0.05, 0.1])
                clf.fit(X, y)
                score = clf.score(X, y)
                # adscore = 1 - (1 - score) * (len(y) - 1) / (len(y) - X.shape[1] - 1)
                # print(len(y), X.shape[1])
                score_list.append(score)
                # print('coef_matrix:\n', clf.coef_)
                # print('model:\n', clf)
                # print('alpha',clf.alpha_ )
                # print('score', score)
                # print('adscore', adscore)

    sum_score = np.sum(score_list)

    ratio_sum = (sumscore_origin - sum_score) / sumscore_origin

    return delete, ratio_sum


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
        for j in range(i + 1, r):
            indices[j] = indices[j - 1] + 1
        yield tuple(pool[i] for i in indices)
