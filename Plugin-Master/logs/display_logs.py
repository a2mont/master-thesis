import matplotlib.pyplot as plt
import numpy as np
import csv
import os
import re

def find_files(regex: re.Pattern, dimension='3D'):
    filenames = []

    for _, _, files in os.walk(f'./{dimension}/'):
        for file in files:
            if regex.match(file):
                filenames.append(file)
    filenames.sort()
    print(f'{len(filenames)} files:\n{filenames}')
    return filenames

def extract_data(files: list, dimension='3D'):
    datas = [[[] for _ in range(len(files))] for _ in range(6)]
    file_id = 0
    for file in files:
        file = f"./{dimension}/" + file
        with open(file) as f:
            reader = csv.reader(f, delimiter=',')
            next(reader)
            for row in reader:
                [datas[i][file_id].append(float(row[i])) for i in range(6) if row[i] != 'nan']
        file_id += 1

    # values for quality is constant
    datas[0] = [datas[0][i][0] for i in range(len(datas[0]))]
    return datas

def plot_quality(q_mins: list, qualities: list):
    rows = int(np.ceil(len(q_mins)/2))
    cols = 2 if len(q_mins) > 1 else 1
    fig, ax = plt.subplots(rows, cols, figsize=(15,15), sharey='all')
    fig.tight_layout(pad=3.5)

    for i,q_min in enumerate(q_mins):
        x = range(0, len(qualities[i]))
        y = qualities[i]
        print(q_min)
        if rows > 1:
            id_row = int(i / 2)
            id_col = i % cols
            ax[id_row, id_col].plot(x,y)
            ax[id_row, id_col].set_title(f'Quality min={q_min}')
            ax[id_row, id_col].text(-1,q_min + 0.01, f'{q_min}', color='r')
            ax[id_row, id_col].axhline(y=float(q_min), color='r', alpha=0.5, linestyle='dashed', label='Minimal quality')
            ax[id_row, id_col].legend()
        elif rows == 1:
            ax[i].plot(x,y)
            ax[i].set_title(f'Quality min={q_min}')
            ax[i].text(-1,q_min + 0.01, f'{q_min}', color='r')
            ax[i].axhline(y=float(q_min), color='r', alpha=0.5, linestyle='dashed', label='Minimal quality')
            ax[i].legend()
        else:
            ax.plot(x,y)
            ax.set_title(f'Quality min={q_min}')
            ax.text(-1,q_min + 0.01, f'{q_min}', color='r')
            ax.axhline(y=float(q_min), color='r', alpha=0.5, linestyle='dashed', label='Minimal quality')
            ax.legend()

    plt.show()

def plot_deltas(q_mins: list, deltas: dict):
    for id,q_min in enumerate(q_mins):
        rows = int(np.ceil(len(deltas)/2))
        cols = 2
        fig, ax = plt.subplots(rows, cols, figsize=(15,15), sharey='all')
        fig.tight_layout(pad=3.5)
        fig.suptitle(f'Quality min = {q_min}')
        n_bins = 23

        for i,(title,delta) in enumerate(deltas.items()):
            data = delta[id]
            id_row = int(i / 2)
            id_col = i % cols
            ax[id_row, id_col].hist(data, ec='darkblue', bins=n_bins)
            ax[id_row, id_col].set_title(title) 

        plt.show()
    
def plot_grouped_deltas(q_mins: list, deltas: dict):
    data_by_quality = {q_min: {title : deltas[title][id] for title in deltas.keys()} for id, q_min in enumerate(q_mins)}
    rows = int(np.ceil(len(q_mins)/2))
    cols = 2 if len(q_mins) > 1 else 1
    fig, ax = plt.subplots(rows, cols, figsize=(15,15), sharey='all')
    fig.tight_layout(pad=3.5)
    n_bins = 23
    for i, q_min in enumerate(q_mins):
        datas = data_by_quality[q_min].values()
        titles = list(data_by_quality[q_min].keys())
        if rows > 1:
            id_row = int(i / 2)
            id_col = i % cols
            ax[id_row, id_col].set_title(f'Quality min={q_min}')
            ax[id_row, id_col].hist(datas, bins=n_bins, histtype='step', stacked=True, fill=False, label=titles)
            ax[id_row, id_col].legend()
        elif rows == 1:
            ax[i].set_title(f'Quality min={q_min}')
            ax[i].hist(datas, bins=n_bins, histtype='step', stacked=True, fill=False, label=titles)
            ax[i].legend()
        else:
            ax.set_title(f'Quality min={q_min}')
            ax.hist(datas, bins=n_bins, histtype='step', stacked=True, fill=False, label=titles)
            ax.legend()
    plt.show()

def main():
    filename = "logs_"
    regex = re.compile(f'({filename}.*\.csv$)')
    files = find_files(regex)
    q_mins,qualities,topos,contras,inserts,smoothes = extract_data(files)
    deltas = {
        "Topological pass": topos,
        "Contraction pass": contras,
        "Insertion pass":   inserts,
        "Smoothing pass":   smoothes}
    plot_quality(q_mins, qualities)
    plot_deltas(q_mins, deltas)
    plot_grouped_deltas(q_mins, deltas)



if __name__ == "__main__":
    main()