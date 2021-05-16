import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv('../MyApplication_h1_Edep.csv', sep=',',header=None, index_col =0)

data.plot(kind='bar')
plt.ylabel('Frequency')
plt.xlabel('Words')
plt.title('Title')

plt.show()
