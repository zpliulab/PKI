import pandas as pd

r=pd.read_csv("results.csv",header=0,index_col=0)
r=r.sort_values(["0"],ascending=False)

lists=[]
for fields in r["1"]:
    fields = fields.strip()
    fields = fields.replace("["," ").replace("]"," ")
    fields = fields.split(",")
    lists.append(fields)
print(lists)

# rel_f=pd.read_csv("rel_f",header=0,index_col=0)
# tf=rel_f["tf"]
# tf=list(set(tf))
# print(len(tf))
# tf=pd.DataFrame(tf)
# tf.to_csv("tf.csv")

all=pd.read_csv("tfuse.csv",index_col=0,header=0)
allcomb=[]
for i in lists:
    comb=[]
    for j in range(0,len(lists[1])):
        #print int(i[j])
        if int(i[j]) == 0:
            print 2
            comb.append(all.iloc[j,0])
    allcomb.append(comb)
print allcomb
allcomb=pd.DataFrame(allcomb)
allcomb.to_csv("tfsort4.csv")

