import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from data_saver.data_loader import DataLoader
from datetime import datetime
import mplcursors

from matplotlib import rc
rc('text', usetex=True)  # Enable LaTeX rendering
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'  # For math formatting

def safe_get(results, key, default="N/A"):
    """Safely get a value from results dictionary with default fallback"""
    return results.get(key, default)

def safe_get_float(results, key, default=0.0): 
    """Safely get a float value with default fallback"""
    try:
        return float(results.get(key, default))
    except (ValueError, TypeError):
        return default

def load_data(filename):
    """Load data from a single file"""
    data_file = os.path.join(os.path.dirname(__file__), 'data', filename)
    
    with DataLoader(data_file) as loader:
        results = loader.load_all()
    
    # Extract all relevant data
    data = {
        'file_name': safe_get(results, 'File_name'),
        'solverType': safe_get(results, 'Alg_unconstrainted'),
        'AlgEqType': safe_get(results, 'Alg_equality'),
        'AlgIneqType': safe_get(results, 'Alg_inequality'),
        'terminationCriterionStr': safe_get(results, 'Termination_criteria:'),
        'linearSolverType': int(safe_get(results, 'linearSolverType', -1)),
        'mpcHorizon': int(safe_get(results, 'mpcHorizon', -1)),
        
        # Signals
        'v_h_signal': np.array(safe_get(results, 'v_h_signal', [])),
        'd_h_signal': np.array(safe_get(results, 'd_h_signal', [])),
        'v_p_signal': np.array(safe_get(results, 'v_p_signal', [])),
        'a_p_signal': np.array(safe_get(results, 'a_p_signal', [])),
        'num_iterations_signal': np.array(safe_get(results, 'num_iterations_signal', [])),
        'cal_time_signal': np.array(safe_get(results, 'cal_time_signal', [])),
        'u_input_signal': np.array(safe_get(results, 'u_input_signal', [])),
        'cost_level': np.array(safe_get(results, 'cost_level', [])),

        
        # AL parameters
        'settingEq0': safe_get_float(results, 'init_lagrange_multiplier_eq', -1),
        'settingEq1': safe_get_float(results, 'init_rho_eq', -1),
        'settingEq2': safe_get_float(results, 'rho_min_eq', -1),
        'settingEq3': safe_get_float(results, 'rho_max_eq', -1),
        'settingEq4': safe_get_float(results, 'rho_update_factor_eq', -1),
        'settingEq5': safe_get_float(results, 'max_num_inner_iterations', -1),
        
        # BIPM parameters
        't_final': safe_get_float(results, 't_final', -1),
        't_init': safe_get_float(results, 't_init', -1),
        'nu_update': safe_get_float(results, 'nu_update', -1)
    }
    
    # Process signals to remove trailing zeros
    if len(data['cal_time_signal']) > 0:
        valid_indices = data['cal_time_signal'] > 0
        data['cal_time_signal'] = data['cal_time_signal'][valid_indices]
        data['num_iterations_signal'] = data['num_iterations_signal'][:len(data['cal_time_signal'])]
        data['v_h_signal'] = data['v_h_signal'][:len(data['cal_time_signal'])]
        data['v_p_signal'] = data['v_p_signal'][:len(data['cal_time_signal'])]
        data['d_h_signal'] = data['d_h_signal'][:len(data['cal_time_signal'])]
        data['u_input_signal'] = data['u_input_signal'][:len(data['cal_time_signal'])]
    
    return data

def main():
    # Set font sizes (IEEE recommends 8-10pt)
   # Set font sizes (IEEE recommends 8-10pt)
    plt.rcParams.update({
        'font.family': 'serif',
        'font.serif': ['Times New Roman'],
        'text.usetex': True,  # Use LaTeX for text rendering
        'font.size': 8,         # Default text size
        'axes.titlesize': 8,    # Title size
        'axes.labelsize': 9,     # Axis labels
        'xtick.labelsize': 8,    # X-axis ticks
        'ytick.labelsize': 8,    # Y-axis ticks
        'legend.fontsize': 8,    # Legend
        'figure.dpi': 600,       # Higher resolution for publications
    })
    
    
    # List of files to load - modify these paths as needed
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
    
    # Load all data
    bipm_data = {horizon: load_data(path) for horizon, path in bipm_files.items()}
    al_data = {horizon: load_data(path) for horizon, path in al_files.items()}
    
    # Create figure with specified size
    scale =.7
    plt.figure(figsize=(4.3 * scale, 2.7 * scale))  # inches
    
    # Plot iterations for each horizon
    horizons = [3, 6, 20]
    bar_width = 0.35
    index = np.arange(len(horizons))
    
    # Calculate means
    bipm_means = [np.mean(bipm_data[h]['num_iterations_signal']) for h in horizons]
    al_means = [np.mean(al_data[h]['num_iterations_signal']) for h in horizons]
    
    # Plot bars
    bipm_bars = plt.bar(index - bar_width/2, bipm_means, bar_width, label='BIPM', color='b')
    al_bars = plt.bar(index + bar_width/2, al_means, bar_width, label='AL', color='r')
    
    # Add value labels on top of each bar
    for bars in [bipm_bars, al_bars]:
        for bar in bars:
            height = bar.get_height()
            plt.text(bar.get_x() + bar.get_width()/2., height,
                     f'{height:.1f}',
                     ha='center', va='bottom')
    
    # Customize plot
    plt.xlabel('MPC Horizon')
    plt.ylabel('Avg. Iterations')
    #plt.title('Comparison of Mean Iterations: BIPM vs AL')
    plt.xticks(index, horizons)
    plt.legend(ncol=2)  # Horizontal legend
    plt.ylim(0, 60)
    #plt.grid(True, axis='y')
    
    plt.tight_layout()
    
    # Save figure
    output_filename = "iterations_bar"
    plt.savefig(f"{output_filename}.pdf", format='pdf', bbox_inches='tight')
    print(f"\nComparison figure saved as: {output_filename}.pdf")
    
    plt.show()

if __name__ == "__main__":
    main()
