import PKIpipeline as PKI
import pandas as pd

networkstructure = pd.read_csv('PKI_submit\\ESCdataset\\Network.csv', header=0)
derivdata = pd.read_csv("PKI_submit\\ESCdataset\\ChangeRate.csv", index_col=0)  # col time
timeseriesdata = pd.read_csv("PKI_submit\\ESCdataset\\TimeSeriesExpressionData.csv", index_col=0)

print(networkstructure.head(), "\n", derivdata.head(), "\n", timeseriesdata.head())
relationship_input, TF = PKI.preprocess2(networkstructure)
sum_score0, mean_socre0 = PKI.Rsquared_origin(relationship=relationship_input, data=timeseriesdata.T, derivdata=derivdata.T)

######################################################################################
resultlist = []
for i in range(0, len(TF)):
    list = [1] * len(TF)
    list[i] = 0
    # print(relationship_input)
    deleteTF, RsquaredSum = PKI.pki(list, TF, relationshipuse=relationship_input, data=timeseriesdata.T,
                                 derivdata = derivdata.T,
                                 sumscore_origin = sum_score0)
    resultlist.append([deleteTF, RsquaredSum])

resultlist = pd.DataFrame(resultlist)
resultlist.columns = ["TF", "pkiscore"]
resultlist = resultlist.sort_values("pkiscore", ascending=False)
print("single_node_rank (for TF in this case):\n",resultlist)
#resultlist.to_csv("resultlist_esc.csv")



