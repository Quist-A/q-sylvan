import os
import argparse
import numpy as np
import matplotlib.pyplot as plt

plt_format = '.png'
bench_path = '../build/benchmark_data/'
combined_file_8 = bench_path + 'combined_8.csv'
group_by = {'grover':'qubits',
            'shor':'qubits',
            'supremacy':'depth',
            'random_circuit':'qubits'}
avg_over = {'grover':'flag',
            'shor':'N',
            'supremacy':'rseed',
            'random_circuit':'rseed'}

args = 0
plot_concur_perf_bool = True
plot_peaknodes_bool   = True
plot_histories_bool   = True
plot_runtimes_bool    = True

def parse_stuff():
    parser = argparse.ArgumentParser()
    parser.add_argument('--replot', help='redo all existing plots',
                        dest='replot', default=False, action='store_true')
    global args 
    args = parser.parse_args()


def plot_history(input_path, output_path, alg_name):
    ys = np.genfromtxt(input_path, dtype=float, delimiter=',', names=True)
    x  = np.arange(start=0, stop=len(ys), step=1)

    _, ax1 = plt.subplots()
    color1 = 'tab:orange'
    color2 = 'tab:blue'

    # gates vs nodes
    ax1.set_xlabel('gates')
    ax1.set_ylabel('qdd nodes', color=color1)
    ax1.plot(x, ys['nodes'], color=color1)
    ax1.tick_params(axis='y', labelcolor=color1)

    # gates vs (total) ctable entries
    ax2 = ax1.twinx()
    ax2.set_ylabel('complex table entries', color=color2)
    ax2.plot(x, ys['amps'], color=color2)
    ax2.tick_params(axis='y', labelcolor=color2)

    plt.title(alg_name.capitalize().replace('_', ' '))
    plt.tight_layout()
    plt.savefig(output_path)
    plt.clf()
    plt.close()


def plot_histories(histories_path, output_folder, alg_name):
    # make sure output folder exists
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # iterate over all .csv files
    for filename in os.listdir(histories_path):
        if filename.endswith('.csv'):
            hist_file_path = histories_path + filename
            output_path = output_folder + filename[:-4] + plt_format
            plot_history(hist_file_path, output_path, alg_name)


def plot_qubits_vs_peak_nodes(data_folder, alg_name):
    input_path  = data_folder + 'summary.csv'
    output_path = data_folder + 'peak_nodes' + plt_format
    data = np.genfromtxt(input_path, dtype=float, delimiter=',', names=True)
    x = data[group_by[alg_name]]
    y = data['peak_nodes']
    plt.scatter(x, y)
    plt.xlabel(group_by[alg_name])
    plt.ylabel('qdd peak nodes')
    plt.title(alg_name.capitalize().replace('_', ' '))
    plt.tight_layout()
    plt.savefig(output_path)
    plt.clf()
    plt.close()

def plot_concurrency_performance(data_folder, alg_name):
    input_path  = data_folder + 'summary.csv'
    output_path = data_folder + 'concurrency' + plt_format
    data = np.genfromtxt(input_path, dtype=float, delimiter=',', names=True)

    lengend_entries = []
    
    # for different number of qubits
    main_group_ids = np.unique(data[group_by[alg_name]])
    #qubits  = np.unique(data['qubits'])
    workers = np.unique(data['workers'])
    for g_id in main_group_ids:
        subset = data[np.where(data[group_by[alg_name]] == g_id)]
        speedups = np.zeros(workers.shape)

        # speedups are calculated withing a group
        group_ids = np.unique(subset[avg_over[alg_name]])
        for group_id in group_ids:
            subsubset = subset[np.where(subset[avg_over[alg_name]] == group_id)]
            speed_w1 = np.mean(subsubset[np.where(subsubset['workers'] == 1)]['runtime'])
            for i, w in enumerate(workers):
                speed_w = np.mean(subsubset[np.where(subsubset['workers'] == w)]['runtime'])
                speedups[i] += (speed_w / speed_w1)**(-1)

        # these are the speedups averaged over the different groups
        speedups /= group_ids.shape

        # actually plot stuff
        plt.plot(workers, speedups)
        w1_times = np.round(subset[np.where(subset['workers'] == 1)]['runtime'], 3)
        leg = '{} {}, time $w_1$ ({},{})'.format(int(g_id), group_by[alg_name], np.min(w1_times), np.max(w1_times))
        lengend_entries.append(leg)

        # add to combined plot if workers are 1,2,4,8
        if (np.all(workers == np.array([1,2,4,8]))):
            with  open(combined_file_8, "a") as f:
                f.write("'{}-{}',".format(alg_name, g_id))
                f.write("{},{},{},{},\n".format(speedups[0], speedups[1], speedups[2], speedups[3]))

    plt.ylabel('average speedup')
    plt.xlabel('number of workers')
    plt.xticks(workers.astype(int))
    plt.legend(lengend_entries)
    plt.plot(workers, np.ones(workers.shape), color='grey', linestyle='--')
    plt.title(alg_name.capitalize().replace('_', ' '))
    plt.savefig(output_path)
    plt.clf()
    plt.close()

def plot_runtimes(data_folder, alg_name, x_axis='qubits'):
    input_path  = data_folder + 'summary.csv'
    output_path = data_folder + 'runtime_vs_' + x_axis + plt_format
    data = np.genfromtxt(input_path, dtype=float, delimiter=',', names=True)

    _, ax1 = plt.subplots()
    color1 = 'tab:orange'
    color2 = 'tab:purple'

    # 'x_axis' vs runtime
    ax1.set_xlabel(x_axis)
    ax1.set_ylabel('runtime', color=color1)
    ax1.scatter(data[x_axis], data['runtime'], color=color1)
    ax1.tick_params(axis='y', labelcolor=color1)

    # 'x_axis' vs avg time per gate
    ax2 = ax1.twinx()
    ax2.set_ylabel('avg time per gate', color=color2)
    ax2.scatter(data[x_axis], data['avg_gate_time'], color=color2)
    ax2.tick_params(axis='y', labelcolor=color2)

    plt.title(alg_name.capitalize().replace('_', ' '))
    plt.tight_layout()
    plt.savefig(output_path)
    plt.clf()
    plt.close()

def init_combined_file():
    with  open(combined_file_8, "w") as f:
        f.write("'name',1,2,4,8\n")


# iterates over all folders in the bench_path and plots everything it can plot
def plot_all():

    if (plot_concur_perf_bool):
        init_combined_file()

    # iterate over all algorithms
    for alg_name in os.listdir(bench_path):

        #alg_path = bench_path + alg_name + '/'
        alg_path = os.path.join(bench_path, alg_name)
        if (not os.path.isdir(alg_path)):
            continue

        # iterate over all experiments
        for exp_folder in os.listdir(alg_path):
            exp_path = os.path.join(alg_path, exp_folder) + '/'
            if (not args.replot and 'concurrency.png' in os.listdir(exp_path)):
                print('skipping {}/{}'.format(alg_name, exp_folder))
            else:
                print('plotting {}/{}'.format(alg_name, exp_folder))

                # plot qubits vs peak nodes
                if (plot_peaknodes_bool):
                    plot_qubits_vs_peak_nodes(exp_path, alg_name)

                # plot concurrency performance
                if (plot_concur_perf_bool):
                    plot_concurrency_performance(exp_path, alg_name)
                
                # plot runtimes vs number of qubits
                if (plot_runtimes):
                    plot_runtimes(exp_path, alg_name, group_by[alg_name])
                    if (alg_name == 'random_circuit'):
                        plot_runtimes(exp_path, alg_name, 'avg_nodes')

                # plot histories for all runs
                if (plot_histories_bool):
                    histories_path = exp_path + 'run_histories/'
                    output_folder  = histories_path + 'plots/'
                    plot_histories(histories_path, output_folder, alg_name)


if __name__ == '__main__':
    parse_stuff()
    plot_all()
