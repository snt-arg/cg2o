import os
import numpy as np
import matplotlib.pyplot as plt
from data_saver.data_loader import DataLoader

from matplotlib import rc
rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

def safe_get(results, key, default="N/A"):
    return results.get(key, default)

def load_data(filename):
    data_file = os.path.join(os.path.dirname(__file__), 'data', filename)
    
    with DataLoader(data_file) as loader:
        results = loader.load_all()
    
    data = {
        'num_iterations_signal': np.array(safe_get(results, 'num_iterations_signal', [])),
        'cal_time_signal': np.array(safe_get(results, 'cal_time_signal', [])),
    }
    
    if len(data['cal_time_signal']) > 0:
        valid_indices = data['cal_time_signal'] > 0
        data['num_iterations_signal'] = data['num_iterations_signal'][valid_indices]
    
    return data

def main():
    plt.rcParams.update({
        'font.size': 10,
        'axes.labelsize': 9,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'legend.fontsize': 8,
        'figure.dpi': 600,
    })
    
    bipm_files = {
        3: "../final/260216_1119_G2O_ISPD_gn_GN_12_3.bin",
        6: "../final/260216_1119_G2O_ISPD_gn_GN_12_3.bin",
        20: "../final/260216_1119_G2O_ISPD_gn_GN_12_3.bin"
    }
    
    al_files = {
        3: "260217_1315_AMPL_0_3.bin",
        6: "260217_1315_AMPL_0_3.bin",
        20: "260217_1315_AMPL_0_3.bin"
    }
    
    bipm_data = {h: load_data(path) for h, path in bipm_files.items()}
    al_data = {h: load_data(path) for h, path in al_files.items()}

    # Create figure with custom grid layout
    fig = plt.figure(figsize=(10, 12))
    gs = fig.add_gridspec(4, 2, height_ratios=[1, 1, 1, 0.7])
    
    # Y-axis limits for each horizon pair
    y_limits = {
        3: (0, 80),
        6: (0, 90),
        20: (0, 300)
    }
    
    # Plot iteration signals
    for i, horizon in enumerate([3, 6, 20]):
        # BIPM plot
        ax1 = fig.add_subplot(gs[i, 0])
        if len(bipm_data[horizon]['num_iterations_signal']) > 0:
            ax1.plot(bipm_data[horizon]['num_iterations_signal'], 'b-', 
                    label=f'BIPM (N={horizon})')
            ax1.set_ylim(y_limits[horizon])
            ax1.set_ylabel('Iterations')
            ax1.legend()
            ax1.grid(True)
        
        # AL plot
        ax2 = fig.add_subplot(gs[i, 1])
        if len(al_data[horizon]['num_iterations_signal']) > 0:
            ax2.plot(al_data[horizon]['num_iterations_signal'], 'r-',
                    label=f'AL (N={horizon})')
            ax2.set_ylim(y_limits[horizon])
            ax2.legend()
            ax2.grid(True)
    
    # Bar plot spanning the entire bottom row
    ax_bar = fig.add_subplot(gs[3, :])
    
    horizons = [3, 6, 20]
    bipm_means = [np.mean(bipm_data[h]['num_iterations_signal']) for h in horizons]
    al_means = [np.mean(al_data[h]['num_iterations_signal']) for h in horizons]
    
    bar_width = 0.35
    index = np.arange(len(horizons))
    
    bipm_bars = ax_bar.bar(index - bar_width/2, bipm_means, bar_width, 
                          label='BIPM', color='b')
    al_bars = ax_bar.bar(index + bar_width/2, al_means, bar_width, 
                        label='AL', color='r')
    
    ax_bar.set_xticks(index)
    ax_bar.set_xticklabels([f'N={h}' for h in horizons])
    ax_bar.set_ylabel('Mean Iterations')
    ax_bar.legend()
    ax_bar.grid(True, axis='y')
    
    # Add value labels on bars
    for bars in [bipm_bars, al_bars]:
        for bar in bars:
            height = bar.get_height()
            ax_bar.text(bar.get_x() + bar.get_width()/2., height,
                       f'{height:.1f}',
                       ha='center', va='bottom', fontsize=8)
    
    plt.tight_layout()
    
    output_filename = "/mnt/e/1.PhD/Papers/paper_BIPM_g2o/figure/MACC/iterations_comparison_grid"
    output_filename = "iterations_comparison_grid"
    plt.savefig(f"{output_filename}.pdf", format='pdf', bbox_inches='tight')
    print(f"Figure saved as: {output_filename}.pdf")
    
    plt.show()

if __name__ == "__main__":
    main()