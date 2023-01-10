import PKIpipeline as PKI
import pandas as pd

networkstructure = pd.read_csv('ESC_submit\\ESCdataset\\Network.csv', header=0)
derivdata = pd.read_csv("ESC_submit\\ESCdataset\\ChangeRate.csv", index_col=0)  # col time
timeseriesdata = pd.read_csv("ESC_submit\\ESCdataset\\TimeSeriesExpressionData.csv", index_col=0)

print(networkstructure.head(), "\n", derivdata.head(), "\n", timeseriesdata.head())
relationship_input, TF = PKI.preprocess2(networkstructure)
sum_score0, mean_socre0 = PKI.Rsquared_origin(relationship=relationship_input, data=timeseriesdata.T, derivdata=derivdata.T)

######################################################################################
two_resultlist=[]
outlist=[]
for i in range(0,len(TF)-1):
    for j in range(i+1,len(TF)):
        outlist.append([i,j])
    print (outlist)


for i in outlist:
    tflist = [1] * len(TF)
    x1=i[0]
    x2=i[1]
    tflist[x1]=0
    tflist[x2]=0
    deleteTF, RsquaredSum = PKI.pki(tflist, TF, relationshipuse=relationship_input, data=timeseriesdata.T, derivdata=derivdata.T,
                               sumscore_origin=sum_score0)
    two_resultlist.append([deleteTF, RsquaredSum])

two_resultlist=pd.DataFrame(two_resultlist)
two_resultlist.columns=["TF","pkiscore"]
two_resultlist=two_resultlist.sort_values("pkiscore",ascending=False)
print("Two TFs combinations:\n",two_resultlist)
#two_resultlist.to_csv("TF2resultlist.csv")


######################################################################################
three_resultlist=[]
for i in list(PKI.combinations(range(len(TF)),3)):
    tflist = [1] * len(TF)
    x1 = i[0]
    x2 = i[1]
    x3 = i[2]
    tflist[x1] = 0
    tflist[x2] = 0
    tflist[x3] = 0
    deleteTF, RsquaredSum = PKI.pki(tflist, TF, relationshipuse=relationship_input, data=timeseriesdata.T, derivdata=derivdata.T,
                               sumscore_origin=sum_score0)
    three_resultlist.append([deleteTF, RsquaredSum])

three_resultlist=pd.DataFrame(three_resultlist)
three_resultlist.columns=["TF","pkiscore"]
three_resultlist=three_resultlist.sort_values(['pkiscore'],ascending=False)
print("Three TFs combinations:\n",three_resultlist)
#three_resultlist.to_csv("TF3resultlist.csv")

#######################################################################################
four_resultlist=[]
for i in list(PKI.combinations(range(len(TF)),4)):
    tflist = [1] * len(TF)
    x1 = i[0]
    x2 = i[1]
    x3 = i[2]
    x4 = i[3]
    tflist[x1] = 0
    tflist[x2] = 0
    tflist[x3] = 0
    tflist[x4] = 0
    deleteTF, RsquaredSum = PKI.pki(tflist, TF, relationshipuse=relationship_input, data=timeseriesdata.T, derivdata=derivdata.T,
                               sumscore_origin=sum_score0)
    four_resultlist.append([deleteTF, RsquaredSum])

four_resultlist=pd.DataFrame(four_resultlist)
four_resultlist.columns=["TF", "pkiscore"]
four_resultlist=four_resultlist.sort_values(['pkiscore'], ascending=False)
print("Four TFs combinations:\n",four_resultlist)
#four_resultlist.to_csv("TF4resultlist.csv")
