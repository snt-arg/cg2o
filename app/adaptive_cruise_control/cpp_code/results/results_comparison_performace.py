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
        'u_input_signal': np.array(safe_get(results, 'cost_level_signal', [])),
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
    
    bipm_filename =  "../final/g2o_13/260217_1506_G2O_ISPD_gn_GN_13_15.bin"
    al_filename =  "../final/IPOPT_10_init/260217_1720_AMPL_0_15.bin"  
    #if len(sys.argv) < 3:
    #    print("Usage: python script.py <BIPM_filename> <AL_filename>")
    #    sys.exit(1)
    

    # Set font sizes (IEEE recommends 8-10pt)
    plt.rcParams.update({
        'font.size': 15,         # Default text size
        'axes.titlesize': 13,    # Title size
        'axes.labelsize': 12,     # Axis labels
        'xtick.labelsize': 12,    # X-axis ticks
        'ytick.labelsize': 12,    # Y-axis ticks
        'legend.fontsize': 12,    # Legend
    })
    
   


 
    # Load data from both files
    bipm_data = load_data(bipm_filename)
    al_data = load_data(al_filename)
   
    bipm_basename = os.path.splitext(os.path.basename(bipm_filename))[0]
    al_basename = os.path.splitext(os.path.basename(al_filename))[0]
    

    
    # Linear solver type mapping
    solver_types = {
        0: "CHOLMOD", 1: "EigenLDLT", 2: "CSparse", 3: "Eigen",
        4: "Dense", 5: "PCG", 6: "EigenLU", 7: "EigenQR",
        8: "EigenSVD", 9: "BiCGSTAB", 10: "EigenDenseQR",
        11: "EigenDenseFullLU", 12: "EigenDensePartialLU",
        13: "EigenPardisoLU", 14: "EigenUmfPackLU", 15: "EigenLDCGlinear"
    }

    # ======================================================================
    # Plotting Section
    # ======================================================================
    num_plots = 5  # Number of subplots
    scale =.7
    plt.figure(figsize=(8.16 * scale, 9 * scale))  # inches
 
    #plt.suptitle('MPC Performance Comparison: BIPM vs AL', y=0.98, fontsize=16)

        # Subplot 1: v_p (blue only)
    ax1 = plt.subplot(num_plots, 1, 1)
    if len(bipm_data['v_p_signal']) > 0:
        ax1.plot(bipm_data['v_p_signal'], 'b-', label=r'$v_p$')  # Using LaTeX math mode
    ax1.set_ylabel(r'$v_p$ (m/s)', fontsize=10)  # Using LaTeX math mode
    ax1.legend()
    ax1.grid(True)
    ax1.set_xticklabels([])


    # Subplot 2: v_h comparison
    ax2 = plt.subplot(num_plots, 1, 2)
    if len(bipm_data['v_h_signal']) > 0:
        ax2.plot(bipm_data['v_h_signal'], 'b-', label='BIPM')
    if len(al_data['v_h_signal']) > 0:
        ax2.plot(al_data['v_h_signal'], 'r--', label='AL')
    ax2.set_ylabel(r'$v_h$ (m/s)')
    ax2.legend()
    ax2.grid(True)
    ax2.set_xticklabels([])

    # Subplot 3: d_h comparison
    ax3 = plt.subplot(num_plots, 1, 3)
    if len(bipm_data['d_h_signal']) > 0:
        ax3.plot(bipm_data['d_h_signal'], 'b-', label='BIPM')
    if len(al_data['d_h_signal']) > 0:
        ax3.plot(al_data['d_h_signal'], 'r--', label='AL')
    ax3.set_ylabel('Distance (m)')
    ax3.legend()
    ax3.grid(True)
    ax3.set_xticklabels([])

    # Subplot 4: u_input_signal comparison
    ax4 = plt.subplot(num_plots, 1, 4)
    if len(bipm_data['u_input_signal']) > 0:
        ax4.plot(bipm_data['u_input_signal'], 'b-', label='BIPM')
    if len(al_data['u_input_signal']) > 0:
        ax4.plot(al_data['u_input_signal'], 'r--', label='AL')
    ax4.set_ylabel('Force Input (N)')
    ax4.set_xlabel('Time Steps')
    ax4.legend()
    ax4.grid(True)

  
    y = al_data['cost_level'] - bipm_data['cost_level']
    # Subplot 5: u_input_signal comparison
    ax5 = plt.subplot(num_plots, 1, 5)
    if len(bipm_data['cost_level']) > 0:
       ax5.plot(y, 'b-', label='BIPM Control')
    if len(al_data['cost_level']) > 0:
       ax5.plot(al_data['cost_level'], 'r--', label='AL Control')
    ax5.set_ylabel('Cost Level')
    ax5.set_xlabel('Time Steps')
    ax5.legend()
    ax5.grid(True)

    plt.tight_layout(rect=[0, 0, 1, 0.96])  # Leave space for suptitle

    # Set font sizes (IEEE recommends 8-10pt)
    plt.rcParams.update({
        'font.family': 'serif',
        'font.serif': ['Times New Roman'],
        'text.usetex': True,  # Use LaTeX for text rendering
        'font.size': 10,         # Default text size
        'axes.titlesize': 10,    # Title size
        'axes.labelsize': 9,     # Axis labels
        'xtick.labelsize': 8,    # X-axis ticks
        'ytick.labelsize': 8,    # Y-axis ticks
        'legend.fontsize': 8,    # Legend
        'figure.dpi': 600,       # Higher resolution for publications
    })
    
    plt.tight_layout(pad=0.5)



    plt.tight_layout(pad=0.5)
    # ======================================================================
    # Text Configuration Section
    # ======================================================================
    if False:  # Set to True to enable text configuration
        text_lines = [
            "Algorithm Configuration Comparison:",
            f"BIPM File: {bipm_data['file_name']} | AL File: {al_data['file_name']}",
            f"Linear Solver: {solver_types.get(bipm_data['linearSolverType'], 'Unknown')} (Type {bipm_data['linearSolverType']}) | Horizon: {bipm_data['mpcHorizon']}",
            "",
            "BIPM Parameters:",
            f"- Inequality: {bipm_data['AlgIneqType']} | Equality: {bipm_data['AlgEqType']} | Unconstrained: {bipm_data['solverType']}",
            f"- t_init: {bipm_data['t_init']:.1f} | t_final: {bipm_data['t_final']:.1e} | ν_update: {bipm_data['nu_update']:.1f}",
            "",
            "AL Parameters:",
            f"- Inequality: {al_data['AlgIneqType']} | Equality: {al_data['AlgEqType']} | Unconstrained: {al_data['solverType']}",
            f"- Eq λ_init: {al_data['settingEq0']:.1e} | ρ_init: {al_data['settingEq1']:.1f} | ρ_min: {al_data['settingEq2']:.1f}",
            f"- ρ_max: {al_data['settingEq3']:.1e} | ρ_update_factor: {al_data['settingEq4']:.1f} | max_inner_iter: {al_data['settingEq5']:.1f}",
            "",
            "Statistics:"
        ]

        # Add iteration statistics
        if len(bipm_data['num_iterations_signal']) > 0 and len(al_data['num_iterations_signal']) > 0:
            text_lines.extend([
                f"- Iterations: BIPM mean={np.mean(bipm_data['num_iterations_signal']):.1f} | AL mean={np.mean(al_data['num_iterations_signal']):.1f}",
                f"- Iterations: BIPM max={np.max(bipm_data['num_iterations_signal'])} | AL max={np.max(al_data['num_iterations_signal'])}"
            ])

        # Add computation time statistics
        if len(bipm_data['cal_time_signal']) > 0 and len(al_data['cal_time_signal']) > 0:
            text_lines.extend([
                f"- Comp Time: BIPM mean={np.mean(bipm_data['cal_time_signal']):.1f}ms | AL mean={np.mean(al_data['cal_time_signal']):.1f}ms",
                f"- Comp Time: BIPM max={np.max(bipm_data['cal_time_signal']):.1f}ms | AL max={np.max(al_data['cal_time_signal']):.1f}ms"
            ])

        plt.figtext(0.5, 0.02, '\n'.join(text_lines), 
                ha='center', fontsize=9, 
                bbox=dict(facecolor='lightgray', alpha=0.5))

        plt.tight_layout()
        plt.subplots_adjust(bottom=0.25)

    # Save figure
    output_filename = "/mnt/e/1.PhD/Papers/paper_BIPM_g2o/figure/MACC"
    output_filename = f"comparison_{bipm_basename}_vs_{al_basename}"

    plt.savefig(f"{output_filename}.pdf", format='pdf', bbox_inches='tight')
    print(f"\nComparison figures saved as: {output_filename}.pdf")
    curser = mplcursors.cursor(hover=False)
    plt.show(block=False)
    plt.pause(5)
    plt.close()

if __name__ == "__main__":
    main()