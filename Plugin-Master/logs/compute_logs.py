import matplotlib.pyplot as plt
import numpy as np
import csv

data = {}
data_smooth = {}
data_no_remesh = {}
idx = 0
with open('quality_logs.csv') as file:
    reader = csv.reader(file, delimiter=',')
    for row in reader:
        if len(row[:-1]) > 0:
            data[idx] = row[:-1]
            idx = idx +1 
with open('quality_logs_smooth_only.csv') as file:
    reader = csv.reader(file, delimiter=',')
    for row in reader:
        if len(row[:-1]) > 0:
            data[idx] = row[:-1]
            idx = idx +1 
with open('quality_logs_nothing.csv') as file:
    reader = csv.reader(file, delimiter=',')
    for row in reader:
        if len(row[:-1]) > 0:
            data_no_remesh[idx] = row[:-1]
            idx = idx +1 

            

data = {id: [float(item) for item in d] for id,d in data.items()}
data_smooth = {id: [float(item) for item in d] for id,d in data_smooth.items()}
data_no_remesh = {id: [float(item) for item in d] for id,d in data_no_remesh.items()}
lengths = [len(d) for d in data.values()]
max_steps = max(lengths)

m = len(data)
n = max(lengths)

print(f'{m} episodes,\nbiggest episode: {n} elements\nlength x: {m*n}')

fig, ax = plt.subplots(figsize=(15,15))

# data = {id: d+[d[-1] for _ in range(len(d),n)] for id,d in data.items()}
x = [(n*i + (j+1)/len(d)) for i,d in data_smooth.items() for j,_ in enumerate(d)]
y = [item for d in data.values() for item in d]
y_smooth = [item for d in data_smooth.values() for item in d]
y_no_remesh = [item for d in data_no_remesh.values() for item in d]

# ax.plot(x,y)  
ax.plot(x,y_smooth)
# ax.plot(x,y_no_remesh)

for i in np.arange(0, m*n, n*2):
    print(i)
    ax.axvspan(i, i + n, facecolor='b', alpha=0.2)
# ax.grid(visible=True)
# ax.set_xticks(np.arange(0, m*n, 1))
# ax.set_yticks(np.arange(0, 1, 0.1))
# ax.axhline(y = 0.65, color='r', linestyle='dashed', label='Quality threshold')
ax.legend(loc='upper right')
    


plt.show()