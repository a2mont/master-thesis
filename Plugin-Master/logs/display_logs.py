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
    datas = names = None
    is_zero = lambda x: abs(float(x)) < 1e-6 
    file_id = 0
    for file in files:
        file = f"./{dimension}/{file}"
        with open(file) as f:
            reader = csv.reader(f, delimiter=',')
            first_line = next(reader)
            names = first_line
            # check for initialization
            if(datas == None):
                datas = [[[] for _ in range(len(files))] for _ in range(len(first_line))]
            for row in reader:
                [datas[i][file_id].append(float(row[i])) if row[i] != 'nan' 
                 else datas[i][file_id].append(0) for i in range(len(first_line)) if i < len(row)]
        file_id += 1
    return datas, names

def inverse_data(data: list):
    return [-float(x) for x in data]


def plot_quality(q_mins: list, qualities: list, deltas: dict, deltas_rejected: dict, title: list = 'Title'):
    rows = int(np.ceil(len(q_mins)/2))
    cols = 2 if len(q_mins) > 1 else 1
    fig, ax = plt.subplots(rows, cols, figsize=(10,8))
    fig.suptitle(title, fontsize=16)
    
    handles = labels = []

    for i,q_min in enumerate(q_mins):
        x = np.arange(0, len(qualities[i])/2, 0.5)
        y = qualities[i]
        if rows > 1:
            id_row = int(i / 2)
            id_col = i % cols
            # twin = ax[id_row,id_col].twinx()
            # plot_deltas(twin, deltas, deltas_rejected)
            quality_ax(ax[id_row,id_col], x, y, -q_min)
            handles,labels = ax[id_row,id_col].get_legend_handles_labels()
        elif cols > 1:
            # twin = ax[i].twinx()
            # plot_deltas(twin, deltas, deltas_rejected)
            quality_ax(ax[i], x, y, -q_min)
            handles,labels = ax[i].get_legend_handles_labels()
        else:
            # twin = ax.twinx()
            # plot_deltas(twin, deltas, deltas_rejected)
            quality_ax(ax, x, y, -q_min)
            handles,labels = ax.get_legend_handles_labels()

    labels_no_dupl = dict(zip(labels, handles))
    fig.legend(labels_no_dupl.values(),labels_no_dupl.keys())


def quality_ax(ax: plt.Axes, x, y ,q_min):
    ax.set_ylim([20,80])
    ax.plot(x,y, label="Final energy")
    ax.set_title(f'Energy threshold = {q_min}')
    ax.text(-1.5,q_min + 0.1, f'{q_min}', color='r')
    ax.axhline(y=float(q_min), color='r', alpha=0.5, linestyle='dashed', label='Energy threshold')

def quality_vector(init: float, q_before: list, q_after: list):
    q_list = [init]
    q_before = [x for x in q_before if x != 0]
    q_after = [x for x in q_after if x != 0]
    pairs = list(zip(q_before,q_after))
    [q_list.extend([before,after]) for before,after in pairs]
    return q_list    

def bar_chart(data, ax: plt.Axes, index = 0, color= '#0000ff'):
    width = 0.25
    offset = width * index
    x = np.arange(len(data))
    colors = color

    copy = np.array([max(min(x, 100000), -100000) for x in data])

    mask_norm = copy > -100000
    mask_inf  = copy <= -100000

    ax.bar(x[mask_norm] + offset,copy[mask_norm], width, color = colors, alpha=0.2)
    # ax.bar(x[mask_inf] + offset ,copy[mask_inf], width, color = colors[1])
    ax.set_yscale('symlog')

def plot_deltas(ax: plt.Axes, deltas: dict, deltas_rejected: dict):
    i = 0
    colors = ['red','blue','green','orange']
    for types in deltas.keys():
        for delta in deltas[types]:
            bar_chart(delta,ax,i,colors[i])
        i += 1

    i = 0
    for types in deltas_rejected.keys():
        for delta in deltas_rejected[types]:
            bar_chart(delta,ax,i,colors[i])
        i += 1


def timesteps_experiment(turn:float=180, q_min: float = 40):
    filename = f"angle_spin/angle{turn}_q_min{q_min}.csv"
    data,names = extract_data([filename])
    angles      = [x for x in data[0][0]]
    times       = [x/1000 for x in data[1][0]]
    quality     = [x for x in data[2][0]]
    quality_avg = [x for x in data[3][0]]
    fig, ax = plt.subplots(2,1, figsize=(10,8), gridspec_kw={'height_ratios': [3, 1]})
    fig.suptitle(f'Angle to time with max deformation = {q_min}')
    ax[0].set_xscale('log')
    ax[1].set_xscale('log')
    ax[0].set_ylabel('Time (s)')
    ax[1].set_xlabel('Number of timesteps')
    ax[1].yaxis.tick_right()
    ax[1].yaxis.set_label_position("right")

    # ax[0].set_ylim([0, 400])

    quality_color = 'red'
    average_color='green'
    base_color = 'blue'

    twin1 = ax[0].twinx()
    # twin1.set_ylim([30,80])
    # twin2 = ax.twinx()
    # ax[1].set_ylim([12.3,12.55])
    twin1.set_ylabel('Deformation energies')

    ax[0].plot(angles,times, color=base_color, label="Execution time")
    twin1.plot(angles, quality, color=quality_color, label="Final worst element")
    twin1.axhline(y=float(q_min), color=quality_color, alpha=0.5, linestyle='dashed', label='Energy threshold')
    ax[1].plot(angles, quality_avg, color=average_color, label="Final mesh average energy")


    lines,labels = ax[0].get_legend_handles_labels()
    lines1,labels1 = twin1.get_legend_handles_labels()
    # lines2,labels2 = ax[1].get_legend_handles_labels()
    ax[0].legend(lines+lines1, labels+labels1)
    ax[1].legend()
    # fig.savefig(f"3D/angle_spin/angle{turn}_q_min{q_min}")

def quality_experiment(experiment):
    filename = "logs_"
    regex = re.compile(f'({filename}.*\.csv$)')
    files = find_files(regex, f"3D/{experiment}")
    if len(files) == 0: return

    data, names = extract_data(files, f"3D/{experiment}")
    data = data[:-1]
    names = names[:-1]
    pairs = zip(names,data)
    label_data = dict((name,value) for name, value in pairs)

    # Q min is constant
    label_data['Quality_min'] = [label_data['Quality_min'][i][0]for i in range(len(files))]
    q_vectors = [quality_vector(label_data['Initial'][i][0], label_data['Before'][i], label_data['After'][i]) for i in range(len(files))]
    q_vectors_avg = [quality_vector(label_data['Initial_avg'][i][0], label_data['Before_avg'][i], label_data['After_avg'][i]) for i in range(len(files))]

    q_vectors = [inverse_data(vec) for vec in q_vectors]
    q_vectors_avg = [inverse_data(vec) for vec in q_vectors_avg]
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
    
    for key,values in deltas.items():
        i = 0
        for value in values:
            deltas[key][i] = [-float(x) if x != 'nan' else 0 for x in value]
            i += 1
    for key,values in deltas_rejected.items():
        i = 0
        for value in values:
            deltas_rejected[key][i] = [-float(x) if x != 'nan' else 0 for x in value]
            i += 1
    
    plot_quality(label_data['Quality_min'], q_vectors, deltas, deltas_rejected, "Deformation of worse element")
    # plot_quality(label_data['Quality_min'], q_vectors_avg, deltas, deltas_rejected, "Average mesh quality")


def turn_experiment(q_min:float = 50):
    filename = f"turn_q_min{q_min}.csv"
    data,names = extract_data([filename], "3D/turn_experiment")
    turn      = [x for x in data[0][0]]
    times       = [x/1000 for x in data[1][0]]
    quality     = [x for x in data[2][0]]
    quality_avg = [x for x in data[3][0]]

    fig, ax = plt.subplots(2,1, figsize=(10,8), gridspec_kw={'height_ratios': [3, 1]})
    fig.suptitle(f'Time to turn with remeshing threshold = {q_min}')

    ax[0].set_xticks(np.arange(90, 450, 90))
    ax[1].set_xticks(np.arange(90, 450, 90))
    ax[0].set_yscale('log')
    ax[0].set_ylabel('Time (s)')
    ax[1].set_xlabel('Total rotation (°)')
    ax[1].yaxis.tick_right()
    ax[1].yaxis.set_label_position("right")

    quality_color = 'red'
    average_color='green'
    base_color = 'blue'

    twin1 = ax[0].twinx()
    twin1.set_ylabel('Deformation energies')

    ax[0].plot(turn,times, color=base_color)
    twin1.plot(turn, quality, color=quality_color, label="Final worst element")
    twin1.axhline(y=float(q_min), color=quality_color, alpha=0.5, linestyle='dashed', label='Remeshing threshold')
    ax[1].plot(turn, quality_avg, color=average_color, label="Final mesh average energy")

    lines,labels = twin1.get_legend_handles_labels()
    lines2,labels2 = ax[1].get_legend_handles_labels()
    ax[0].legend(lines, labels, loc='upper left')
    ax[1].legend(lines2, labels2, loc='upper left')

def wm_experiment():
    filename = "logs_"
    regex = re.compile(f'({filename}.*\.csv$)')
    files = find_files(regex, f"3D/wm_experiment")
    if len(files) == 0: return

    data, names = extract_data(files, f"3D/wm_experiment")
    data = data[:-1]
    names = names[:-1]
    pairs = zip(names,data)
    label_data = dict((name,value) for name, value in pairs)

    # Q min is constant
    label_data['Quality_min'] = [-label_data['Quality_min'][i][0] for i in range(len(files))]
    q_vectors = [quality_vector(label_data['Initial'][i][0], label_data['Before'][i], label_data['After'][i]) for i in range(len(files))]
    q_vectors = [inverse_data(vec) for vec in q_vectors]

    fig,ax = plt.subplots(2, figsize=(10,8))
    x = np.arange(0, len(q_vectors[0])/2, 0.5)
    ax[0].plot(x,q_vectors[0], label="No world mesh")
    ax[0].plot(x,q_vectors[1], label="With world mesh")
    # ax[0].text(-1.5,label_data['Quality_min'][0] + 0.1, f'{label_data['Quality_min'][0]}', color='r')
    ax[0].axhline(y=float(label_data['Quality_min'][0]), color='r', alpha=0.5, linestyle='dashed', label='Max energy')
    
    ax[1].plot(x,q_vectors[0])
    ax[1].plot(x,q_vectors[1])
    # ax[1].text(-1.5,label_data['Quality_min'][0] + 0.1, f'{label_data['Quality_min'][0]}', color='r')
    ax[1].axhline(y=float(label_data['Quality_min'][0]), color='r', alpha=0.5, linestyle='dashed')
    ax[1].set_xlim(180,360)
    fig.legend()

def stretch_experiment(stretch_coeff, q_min):    
    filename = f"stretch_experiment/stretch{stretch_coeff}_qmin_{q_min}.csv"
    data,names = extract_data([filename])
    angles      = [x for x in data[0][0]]
    times       = [x/1000 for x in data[1][0]]
    quality     = [x for x in data[2][0]]
    quality_avg = [x for x in data[3][0]]

    fig, ax = plt.subplots(2,1, figsize=(10,8), gridspec_kw={'height_ratios': [3, 1]})
    fig.suptitle(f'Stretch to time with remeshing threshold = {q_min}')
    ax[0].set_xscale('log')
    ax[1].set_xscale('log')
    ax[0].set_ylabel('Time (s)')
    ax[1].set_xlabel('Timesteps')

    ax[0].set_ylim([0, 120])

    quality_color = 'red'
    average_color='green'
    base_color = 'blue'

    twin1 = ax[0].twinx()
    twin1.set_ylim([30, 58])
    # twin2 = ax.twinx()
    ax[1].set_ylim([12.5, 14])
    ax[1].set_ylabel('Deformation energies')
    ax[1].yaxis.tick_right()
    ax[1].yaxis.set_label_position("right")

    ax[0].plot(angles,times, color=base_color, label='Time')
    twin1.plot(angles, quality, color=quality_color, label="Worst element")
    twin1.axhline(y=float(q_min), color=quality_color, alpha=0.5, linestyle='dashed', label='Maximal energy')
    ax[1].plot(angles, quality_avg, color=average_color, label="Mesh average energy")

    lines,labels = ax[0].get_legend_handles_labels()
    lines2,labels2 = twin1.get_legend_handles_labels()
    ax[0].legend(lines+lines2, labels+labels2)
    ax[1].legend()
    fig.savefig(f"3D/stretch_experiment/stretch{q_min}")
    

def length_experiment(q_min:float = 50):
    filename = f"length_q_min{q_min}.csv"
    data,names = extract_data([filename], "3D/length_experiment")
    length      = [x for x in data[0][0]]
    times       = [x/1000 for x in data[1][0]]
    quality     = [x for x in data[2][0]]
    quality_avg = [x for x in data[3][0]]

    fig, ax = plt.subplots(2,1, figsize=(10,8), gridspec_kw={'height_ratios': [3, 1]})
    fig.suptitle(f'Time to length with remeshing threshold = {q_min}')

    ax[0].set_xticks(length)
    ax[1].set_xticks(length)
    ax[0].set_ylabel('Time (s)')
    ax[1].set_xlabel('Stretching factor')

    quality_color = 'red'
    average_color='green'
    base_color = 'blue'

    twin1 = ax[0].twinx()
    twin1.tick_params(axis='y', colors=quality_color)
    # twin2 = ax.twinx()
    ax[1].tick_params(axis='y', colors=average_color)
    ax[1].set_ylabel('Deformation energies')
    ax[1].yaxis.tick_right()
    ax[1].yaxis.set_label_position("right")

    ax[0].plot(length,times, color=base_color)
    twin1.plot(length, quality, color=quality_color, label="Final worst element")
    twin1.axhline(y=float(q_min), color=quality_color, alpha=0.5, linestyle='dashed', label='Remeshing threshold')
    ax[1].plot(length, quality_avg, color=average_color, label="Final mesh average energy")
    
    lines,labels = ax[0].get_legend_handles_labels()
    lines2,labels2 = twin1.get_legend_handles_labels()
    ax[0].legend(lines+lines2, labels+labels2, loc='upper left')
    ax[1].legend()
    fig.savefig(f"3D/length_experiment/length{q_min}")

    


def main():
    # Angle timesteps experiment
    # timesteps_experiment(180,35)
    # timesteps_experiment(180,40)
    # timesteps_experiment(180,50)
    # timesteps_experiment(360,80)

    # Quality spin epxeriment
    experiments = ["quality_spin","aquality_stretch","arevert_experiment"]
    # [quality_experiment(exp) for exp in experiments]
    

    # Turn experiment
    # turn_experiment()

    # Length experiment
    # length_experiment(40)

    # Stretch timesteps experiment
    # stretch_experiment("1_7x",35)
    # stretch_experiment("1_7x",40)

    # Wm vs no wm experiment
    # wm_experiment()

    plt.show()
    



if __name__ == "__main__":
    main()



# ------------------- OLD ------------------------------------

# def plot_deltas(q_mins: list, deltas: dict):
#     for id,q_min in enumerate(q_mins):
#         rows = int(np.ceil(len(deltas)/2))
#         cols = 2
#         fig, ax = plt.subplots(rows, cols, figsize=(15,15), sharey='all', sharex='all')
#         fig.tight_layout(pad=3.5)
#         fig.suptitle(f'Quality min = {q_min}')
#         n_bins = 23

#         for i,(title,delta) in enumerate(deltas.items()):
#             data = delta[id]
#             id_row = int(i / 2)
#             id_col = i % cols
#             ax[id_row, id_col].hist(data, ec='darkblue', bins=n_bins)
#             ax[id_row, id_col].set_title(title) 

#         plt.show()
    
# def plot_grouped_deltas(q_mins: list, deltas: dict):
#     data_by_quality = {q_min: {title : deltas[title][id] for title in deltas.keys()} for id, q_min in enumerate(q_mins)}
#     rows = int(np.ceil(len(q_mins)/2))
#     cols = 2 if len(q_mins) > 1 else 1
#     fig, ax = plt.subplots(rows, cols, figsize=(15,15), sharey='all', sharex='all')
#     fig.tight_layout(pad=3.5)
#     n_bins = 23
#     for i, q_min in enumerate(q_mins):
#         datas = data_by_quality[q_min].values()
#         print(datas)
#         titles = list(data_by_quality[q_min].keys())
#         if rows > 1:
#             id_row = int(i / 2)
#             id_col = i % cols
#             ax[id_row, id_col].set_title(f'Quality min={q_min}')
#             ax[id_row, id_col].hist(datas, bins=n_bins, histtype='step', stacked=True, fill=False, label=titles)
#             ax[id_row, id_col].legend()
#         elif cols > 1:
#             ax[i].set_title(f'Quality min={q_min}')
#             for data in datas:
#                 x = np.arange(len(data))
#                 ax[i].bar(x,data, label='AAAA')
#             ax[i].legend()
#         else:
#             ax.set_title(f'Quality min={q_min}')
#             ax.hist(datas, bins=n_bins, histtype='step', stacked=True, fill=False, label=titles)
#             ax.legend()
#     plt.show()

