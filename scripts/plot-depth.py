import matplotlib.pyplot as plt
import pandas as pd


in_file = snakemake.input[0]
out_file = snakemake.output[0]

df_file = pd.read_csv(in_file, sep='\t', header=None, names=['col1','col2','col3'])  


df1 = df_file[df_file['col1'] == "chr1H_LR890096.1"]
df2 = df_file[df_file['col1'] == "chr2H_LR890097.1"]
df3 = df_file[df_file['col1'] == "chr3H_LR890098.1"]
df4 = df_file[df_file['col1'] == "chr4H_LR890099.1"]
df5 = df_file[df_file['col1'] == "chr5H_LR890100.1"]
df6 = df_file[df_file['col1'] == "chr6H_LR890101.1"]
df7 = df_file[df_file['col1'] == "chr7H_LR890102.1"]
df8 = df_file.loc[373185:382184]


fig, axs = plt.subplots(4,2,sharex=True, sharey=True)

axs[0,0].plot("col2",
         "col3",
         'tab:orange',
         data=df1)

axs[0,1].plot("col2",
         "col3",
         'tab:green',
         data=df2)

axs[1,0].plot("col2",
         "col3",'tab:red',
         data=df3)     

axs[1,1].plot("col2",
         "col3",'tab:blue',
         data=df4)

axs[2,0].plot("col2",
         "col3",'tab:brown',
         data=df5)

axs[2,1].plot("col2",
         "col3",'tab:pink',
         data=df6)

axs[3,0].plot("col2",
         "col3",'tab:grey',
         data=df7)
axs[3,1].plot("col2",
         "col3",'tab:cyan',
         data=df8)

print (df8.shape)
plt.savefig(out_file)
