# import matplotlib.pyplot as plt
# from pandas.tools.plotting import parallel_coordinates
# from pandas import read_csv
# 
# sample = read_csv('iris.data.txt')
# plt.figure()
# parallel_coordinates(sample, 'Name')

from pandas import read_csv
from pandas.tools.plotting import parallel_coordinates
from matplotlib import pyplot as plt
df = read_csv('iris.data')
parallel_coordinates(df, 'Name', color=('#556270','#4ECDC4', '#C7F464'))
plt.show()