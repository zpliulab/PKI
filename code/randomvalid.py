import random
import pandas as pd

tf2=pd.read_csv("tfsort2new.csv",header=0,index_col=0)
tf3=pd.read_csv("tfsort3new.csv",header=0,index_col=0)
randomtf2=[random.randint(0,len(tf2)) for _ in range(100)] #tf2 Index
randomtf3=[random.randint(0,len(tf3)) for _ in range(100)]
tf2best= "('LEF1', 'SOX2')"
for i in randomtf2:
    tf2.iloc




