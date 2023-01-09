import PKIpipeline as PKI
import pandas as pd

# read file
networkstructure = pd.read_csv('PKI_submit\DreamDataset\DREAM3GoldStandard_InSilicoSize100_Ecoli1.txt',
                               sep="\t", header=None)
derivdata = pd.read_csv("PKI_submit\DreamDataset\ChangeRate21_dream3E1_100.csv")  # col time
timeseriesdata = pd.read_csv("PKI_submit\DreamDataset\InSilicoSize100-Ecoli1-trajectories.tsv",
                             sep="\t", index_col=0)


relationship_input, data_input, derivdata_input, TF = PKI.preprocess(networkstructure, timeseriesdata, derivdata)
sum_score_based_changeRate, mean_score_based_changeRate = PKI.Rsquared_origin(relationship=relationship_input,
                                                                       data=data_input,
                                                                       derivdata=derivdata_input)
# compare with diff
da1 = data_input.iloc[1:21, ]
da2 = data_input.iloc[0:20, ]
da2.index = da1.index
diffdata = da1 - da2

sum_score_based_diff, mean_score_based_diff = PKI.Rsquared_origin(relationship=relationship_input,
                                                           data=data_input.iloc[1:21, ],
                                                           derivdata=diffdata)

print("*R-squared value*","\n-based on ODE model in PKI (Sum, Mean): ", sum_score_based_changeRate, mean_score_based_changeRate,
      "\n-based on difference method model (Sum, Mean): ", sum_score_based_diff, mean_score_based_diff)
