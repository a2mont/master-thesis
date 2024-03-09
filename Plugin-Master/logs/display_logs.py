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

def extract_data(files: list, drop_zeros: bool = False, dimension='3D'):
    datas = None
    is_zero = lambda x: abs(float(x)) < 1e-6 
    file_id = 0
    for file in files:
        file = f"./{dimension}/{file}"
        with open(file) as f:
            reader = csv.reader(f, delimiter=',')
            first_line = next(reader)
            # check for initialization
            if(datas == None):
                datas = [[[] for _ in range(len(files))] for _ in range(len(first_line[:-1]))]
            for row in reader:
                if drop_zeros:
                    [datas[i][file_id].append(float(row[i])) for i in range(len(first_line[:-1])) if row[i] != 'nan' and not is_zero(row[i])]
                else:
                    [datas[i][file_id].append(float(row[i])) for i in range(len(first_line[:-1])) if row[i] != 'nan']
        file_id += 1

    return datas, first_line[:-1]

def treat_data(input: list, tranform_nan: bool):
    output = []
    if(tranform_nan):
        output = [x if x != 'nan' else 0 for x in input]
    else:
        output = [x for x in input if x != 'nan']
    return output


def plot_quality(q_mins: list, qualities: list, baseline: list, title: list = 'Title'):
    rows = int(np.ceil(len(q_mins)/2))
    cols = 2 if len(q_mins) > 1 else 1
    fig, ax = plt.subplots(rows, cols, figsize=(15,15), sharey='all')
    plt.xticks(np.arange(0, len(baseline), step=1))
    # fig.tight_layout(pad=3.5)
    fig.suptitle(title, fontsize=16)

    for i,q_min in enumerate(q_mins):
        if len(qualities[i]) != len(baseline):
            print(f'Quality {q_min} size does not match baseline size: {len(qualities[i])} vs {len(baseline)}')
            continue
        x = np.arange(0, len(qualities[i])/2, 0.5)
        y = qualities[i]
        print(y)
        print(baseline)
        if rows > 1:
            id_row = int(i / 2)
            id_col = i % cols
            quality_ax(ax[id_row,id_col], x, y, baseline, q_min)
        elif cols > 1:
            quality_ax(ax[i], x, y, baseline, q_min)
        else:
            quality_ax(ax, x, y, baseline, q_min)

    plt.show()

def quality_ax(ax: plt.Axes, x, y, baseline ,q_min):
    step = 0.5
    ax.plot(x,y, label='With remeshing')
    ax.plot(x,baseline, label='No remeshing',color='g')
    for i in np.arange(0, len(baseline)/2, 2*step):
        ax.axvspan(i, i+step, facecolor='r', alpha=0.2)
    for i in np.arange(0, len(baseline)/2, 2*step):
        ax.axvspan(i+step, i+step*2, facecolor='b', alpha=0.2)
    ax.set_title(f'Quality min={q_min}')
    ax.text(-1.5,q_min + 0.1, f'{q_min}', color='r')
    ax.axhline(y=float(q_min), color='r', alpha=0.5, linestyle='dashed', label='Minimal quality')
    ax.legend()

def plot_deltas(q_mins: list, deltas: dict):
    for id,q_min in enumerate(q_mins):
        rows = int(np.ceil(len(deltas)/2))
        cols = 2
        fig, ax = plt.subplots(rows, cols, figsize=(15,15), sharey='all', sharex='all')
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
    
# def plot_grouped_deltas(q_mins: list, deltas: dict):
#     data_by_quality = {q_min: {title : deltas[title][id] for title in deltas.keys()} for id, q_min in enumerate(q_mins)}
#     rows = int(np.ceil(len(q_mins)/2))
#     cols = 2 if len(q_mins) > 1 else 1
#     fig, ax = plt.subplots(rows, cols, figsize=(15,15), sharey='all', sharex='all')
#     fig.tight_layout(pad=3.5)
#     n_bins = 23
#     for i, q_min in enumerate(q_mins):
#         datas = data_by_quality[q_min].values()
#         titles = list(data_by_quality[q_min].keys())
#         if rows > 1:
#             id_row = int(i / 2)
#             id_col = i % cols
#             ax[id_row, id_col].set_title(f'Quality min={q_min}')
#             ax[id_row, id_col].hist(datas, bins=n_bins, histtype='step', stacked=True, fill=False, label=titles)
#             ax[id_row, id_col].legend()
#         elif cols > 1:
#             ax[i].set_title(f'Quality min={q_min}')
#             ax[i].hist(datas, bins=n_bins, histtype='step', stacked=True, fill=False, label=titles)
#             ax[i].legend()
#         else:
#             ax.set_title(f'Quality min={q_min}')
#             ax.hist(datas, bins=n_bins, histtype='step', stacked=True, fill=False, label=titles)
#             ax.legend()
#     plt.show()

def plot_grouped_deltas(q_mins: list, deltas: dict):
    data_by_quality = {q_min: {title : deltas[title][id] for title in deltas.keys()} for id, q_min in enumerate(q_mins)}
    rows = int(np.ceil(len(q_mins)/2))
    cols = 2 if len(q_mins) > 1 else 1
    fig, ax = plt.subplots(rows, cols, figsize=(15,15), sharey='all', sharex='all')
    fig.tight_layout(pad=3.5)
    n_bins = 23
    for i, q_min in enumerate(q_mins):
        datas = data_by_quality[q_min].values()
        print(datas)
        titles = list(data_by_quality[q_min].keys())
        if rows > 1:
            id_row = int(i / 2)
            id_col = i % cols
            ax[id_row, id_col].set_title(f'Quality min={q_min}')
            ax[id_row, id_col].hist(datas, bins=n_bins, histtype='step', stacked=True, fill=False, label=titles)
            ax[id_row, id_col].legend()
        elif cols > 1:
            ax[i].set_title(f'Quality min={q_min}')
            for data in datas:
                x = np.arange(len(data))
                ax[i].bar(x,data, label='AAAA')
            ax[i].legend()
        else:
            ax.set_title(f'Quality min={q_min}')
            ax.hist(datas, bins=n_bins, histtype='step', stacked=True, fill=False, label=titles)
            ax.legend()
    plt.show()


def quality_vector(init: float, q_before: list, q_after: list):
    q_list = [init]
    pairs = list(zip(q_before,q_after))
    [q_list.extend([before,after]) for before,after in pairs]
    return q_list    

def bar_chart(data, ax: plt.Axes, index = 0, color= '#0000ff'):
    width = 0.25
    offset = width * index
    x = np.arange(len(data))
    colors = (color, color + '50')

    copy = np.array([y if y > -1000 else -100000 for y in data])

    mask_norm = copy > -100000
    mask_inf  = copy <= -100000

    ax.bar(x[mask_norm] + offset,copy[mask_norm], width, color = colors[0])
    ax.bar(x[mask_inf] + offset ,copy[mask_inf], width, color = colors[1])
    ax.set_yscale('symlog')

def main():
    filename = "logs_"
    regex = re.compile(f'({filename}.*\.csv$)')
    files = find_files(regex)
    data, names = extract_data(files)
    pairs = zip(names,data)
    label_data = dict((name,value) for name, value in pairs)

    baseline_data, baseline_names = extract_data(['baseline_20x20_20t.csv'])
    baseline_pairs = zip(baseline_names,baseline_data)
    baseline = dict((name,value) for name, value in baseline_pairs)

    baseline_quality = quality_vector(baseline['Initial'][0][0], baseline['Before'][0], baseline['After'][0])
    baseline_quality_average = quality_vector(baseline['Initial_avg'][0][0], baseline['Before_avg'][0], baseline['After_avg'][0])

    # Q min is constant
    label_data['Quality_min'] = [label_data['Quality_min'][i][0]for i in range(len(files))]

    q_vectors = [quality_vector(label_data['Initial'][i][0], label_data['Before'][i], label_data['After'][i]) for i in range(len(files))]
    q_vectors_avg = [quality_vector(label_data['Initial_avg'][i][0], label_data['Before_avg'][i], label_data['After_avg'][i]) for i in range(len(files))]
    deltas = {
        "Topological pass": label_data['Topological'],
        "Contraction pass": label_data['Contraction'],
        "Insertion pass":   label_data['Insertion'],
        "Smoothing pass":   label_data['Smoothing']}
    deltas_rejected = {
        "Topological pass": label_data['Topological_reject'],
        "Contraction pass": label_data['Contraction_reject'],
        "Insertion pass":   label_data['Insertion_reject'],
        "Smoothing pass":   label_data['Smoothing_reject']}
    plot_quality(label_data['Quality_min'], q_vectors, baseline_quality, "Quality of worse element")
    plot_quality(label_data['Quality_min'], q_vectors_avg, baseline_quality_average, "Average mesh quality")
    # plot_deltas(q_mins, deltas)
    # plot_grouped_deltas(label_data['Quality_min'], deltas)

    fig, ax = plt.subplots(1, 1, figsize=(12,10))
    bar_chart(deltas['Topological pass'][0], ax, 0, '#003f5c')
    bar_chart(deltas['Contraction pass'][0], ax, 1, '#bc5090')
    bar_chart(deltas['Insertion pass'][0], ax, 2, '#ffa600')
    bar_chart(deltas['Smoothing pass'][0], ax, 3, '#ffa600')
    bar_chart(deltas_rejected['Topological pass'][0], ax, 0, '#003f5c')
    bar_chart(deltas_rejected['Contraction pass'][0], ax, 1, '#bc5090')
    bar_chart(deltas_rejected['Insertion pass'][0], ax, 2, '#ffa600')
    bar_chart(deltas_rejected['Smoothing pass'][0], ax, 2, '#ffa600')
    plt.show()



if __name__ == "__main__":
    main()