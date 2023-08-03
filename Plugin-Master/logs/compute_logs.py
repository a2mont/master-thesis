import matplotlib.pyplot as plt
import numpy as np
import csv
import os
import re


def plot_logs(datas):
    rows = int(np.ceil(len(datas)/2))
    cols = 2
    fig, ax = plt.subplots(rows, cols, figsize=(15,15), sharey='all')
    fig.tight_layout(pad=3.5)

    for i, (title, model) in enumerate(datas.items()):
        if len(model) == 0:
            continue
        id_row = int(i / 2)
        id_col = i % cols
        lengths = [len(d) for d in model.values()]
        m = len(model)
        n = max(lengths) 
        x = [(n*id + (j+1)/len(dat)) for id,dat in model.items() for j,_ in enumerate(dat)]
        x = [val - min(x) for val in x]
        y = [item for d in model.values() for item in d]
        if rows > 1:
            ax[id_row, id_col].plot(x,y)
            ax[id_row, id_col].set_title(f'Quality min={title}')
            ax[id_row, id_col].axhline(y=float(title), color='r', alpha=0.5, linestyle='dashed', label='Q_min')
            ax[id_row, id_col].legend()
            for j in np.arange(0, m*n, n*2):
                ax[id_row, id_col].axvspan(j, j + n, facecolor='b', alpha=0.2)
            # if not ax[-1,1].lines: ax[-1,1].set_visible(False)
        else:
            ax[i].plot(x,y)
            ax[i].set_title(f'Quality min={title}')
            ax[i].axhline(y=float(title), color='r', alpha=0.5, linestyle='dashed', label='Q_min')
            ax[i].legend()
            for j in np.arange(0, m*n, n*2):
                ax[i].axvspan(j, j + n, facecolor='b', alpha=0.2)
            # if not ax[-1].lines: ax[-1].set_visible(False)
        # print(f'{i}:\nm = {m},\nn = {n} \nlength x: {m*n}')
    plt.show()

def plot_timesteps(datas):
    datas = {key: data[1] for key, data in datas.items() }
    rows = int(np.ceil(len(datas)/2))
    cols = 2
    fig, ax = plt.subplots(rows, cols, figsize=(15,15), sharey='all')
    fig.tight_layout(pad=3.5)
    for i, (title, model) in enumerate(datas.items()):
        if len(model) == 0:
            continue
        id_row = int(i / 2)
        id_col = i % cols
        x = range(0, len(model))
        y = model
        if rows > 1:
            ax[id_row, id_col].plot(x,y)
            ax[id_row, id_col].set_title(f'Quality min={title}')
            ax[id_row, id_col].axhline(y=float(title), color='r', alpha=0.5, linestyle='dashed', label='Q_min')
            ax[id_row, id_col].legend()
            # if not ax[-1,1].lines: ax[-1,1].set_visible(False)
        else:
            ax[i].plot(x,y)
            ax[i].set_title(f'Quality min={title}')
            ax[i].axhline(y=float(title), color='r', alpha=0.5, linestyle='dashed', label='Q_min')
            ax[i].legend()
            # if not ax[-1].lines: ax[-1].set_visible(False)

    plt.show()

def find_files(regex):
    filenames = []

    for _, _, files in os.walk('.'):
        for file in files:
            if regex.match(file):
                filenames.append(file)
    filenames.sort()
    print(f'{len(filenames)} files\n{filenames}')
    return filenames

def data_from_files(filenames):
    datas = {}
    for f in filenames:
        idx = 0
        with open(f) as file:
            reader = csv.reader(file, delimiter=',')
            name = ''
            for row in reader:
                if reader.line_num == 1:
                    continue
                if reader.line_num == 2:
                    name = row[0] 
                    datas[name] = {}
                datas[name][idx] = [float(item) for item in row[:-1]]
                idx = idx + 1
    return datas 

def main():
    # Constant logs
    regex = re.compile('(.*logs.*\.csv$)')
    filenames = find_files(regex)
    datas_logs = data_from_files(filenames)
    plot_logs(datas_logs)

    # Timestep logs
    regex = re.compile('(.*timesteps.*\.csv$)')
    filenames = find_files(regex)
    datas_ts = data_from_files(filenames)
    plot_timesteps(datas_ts)




if __name__ == "__main__":
    main()