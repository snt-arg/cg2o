import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from data_saver.data_loader import DataLoader
from datetime import datetime
import mplcursors

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
    #if len(sys.argv) < 3:
    #    print("Usage: python script.py <BIPM_filename> <AL_filename>")
    #    sys.exit(1)
    
    bipm_filename =  "accepted/250526_1844_BIPM_gn_GN_13_6.bin"
    al_filename = "accepted/250528_1302_AL_al_GN_13_6.bin"
   
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
    plt.figure(figsize=(15, 12))
    # Enable interactive mode
    plt.ion()
    plt.suptitle('MPC Performance Comparison: BIPM vs AL', y=0.98, fontsize=16)

    # Subplot 1: v_p and v_h (BIPM vs AL)
    ax1 = plt.subplot(3, 2, 1)
    if len(bipm_data['v_p_signal']) > 0:
        ax1.plot(bipm_data['v_p_signal'], 'b-', label='v_p')
    if len(bipm_data['v_h_signal']) > 0:
        ax1.plot(bipm_data['v_h_signal'], 'b-', label='BIPM: v_h')
    if len(al_data['v_h_signal']) > 0:
        ax1.plot(al_data['v_h_signal'], 'r-', label='AL: v_h')
    ax1.set_ylabel('Velocity (m/s)')
    ax1.legend()
    ax1.grid(True)

    # Subplot 2: Iterations comparison
    ax2 = plt.subplot(3, 2, 2)
    if len(bipm_data['num_iterations_signal']) > 0:
        ax2.plot(bipm_data['num_iterations_signal'], 'b-', label='BIPM Iterations')
        bipm_mean = np.mean(bipm_data['num_iterations_signal'])
        ax2.axhline(bipm_mean, color='b', linestyle='--', label=f'BIPM Mean: {bipm_mean:.1f}')
    if len(al_data['num_iterations_signal']) > 0:
        ax2.plot(al_data['num_iterations_signal'], 'r-', label='AL Iterations')
        al_mean = np.mean(al_data['num_iterations_signal'])
        ax2.axhline(al_mean, color='r', linestyle='--', label=f'AL Mean: {al_mean:.1f}')
    ax2.set_ylabel('Iterations')
    ax2.legend()
    ax2.grid(True)

    # Subplot 3: Distance comparison
    ax3 = plt.subplot(3, 2, 3)
    if len(bipm_data['d_h_signal']) > 0:
        ax3.plot(bipm_data['d_h_signal'], 'b-', label='BIPM Distance')
    if len(al_data['d_h_signal']) > 0:
        ax3.plot(al_data['d_h_signal'], 'r-', label='AL Distance')
    ax3.set_ylabel('Inter-distance (m)')
    ax3.legend()
    ax3.grid(True)

    # Subplot 4: Computation time comparison
    ax4 = plt.subplot(3, 2, 4)
    if len(bipm_data['cal_time_signal']) > 0:
        ax4.plot(bipm_data['cal_time_signal'], 'b-', label='BIPM Time')
        bipm_time_mean = np.mean(bipm_data['cal_time_signal'])
        ax4.axhline(bipm_time_mean, color='b', linestyle='--', label=f'BIPM Mean: {bipm_time_mean:.1f}ms')
    if len(al_data['cal_time_signal']) > 0:
        ax4.plot(al_data['cal_time_signal'], 'r-', label='AL Time')
        al_time_mean = np.mean(al_data['cal_time_signal'])
        ax4.axhline(al_time_mean, color='r', linestyle='--', label=f'AL Mean: {al_time_mean:.1f}ms')
    ax4.set_ylabel('Time (ms)')
    ax4.legend()
    ax4.grid(True)

    # Subplot 5: Control input comparison
    ax5 = plt.subplot(3, 2, 5)
    if len(bipm_data['u_input_signal']) > 0:
        ax5.plot(bipm_data['u_input_signal'], 'b-', label='BIPM Control')
    if len(al_data['u_input_signal']) > 0:
        ax5.plot(al_data['u_input_signal'], 'r-', label='AL Control')
    ax5.set_xlabel('Time steps')
    ax5.set_ylabel('Control Input (N)')
    ax5.legend()
    ax5.grid(True)

    # Subplot 6: Time per iteration comparison

    ax6 = plt.subplot(3, 2, 6)
    if len(bipm_data['cal_time_signal']) > 0 and len(bipm_data['num_iterations_signal']) > 0:
        bipm_time_per_iter = bipm_data['cal_time_signal'] / bipm_data['num_iterations_signal']
        ax6.plot(bipm_time_per_iter, 'b-', label='BIPM Time/Iter')
        bipm_tpi_mean = np.mean(bipm_time_per_iter)
        ax6.axhline(bipm_tpi_mean, color='b', linestyle='--', label=f'BIPM Mean: {bipm_tpi_mean:.2f}ms')
    if len(al_data['cal_time_signal']) > 0 and len(al_data['num_iterations_signal']) > 0:
        al_time_per_iter = al_data['cal_time_signal'] / al_data['num_iterations_signal']
        ax6.plot(al_time_per_iter, 'r-', label='AL Time/Iter')
        al_tpi_mean = np.mean(al_time_per_iter)
        ax6.axhline(al_tpi_mean, color='r', linestyle='--', label=f'AL Mean: {al_tpi_mean:.2f}ms')
    ax6.set_xlabel('Time steps')
    ax6.set_ylabel('Time/Iter (ms)')
    ax6.legend()
    ax6.grid(True)

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
    output_filename = f"comparison_{bipm_basename}_vs_{al_basename}"
    plt.savefig(f"{output_filename}.pdf", format='pdf', bbox_inches='tight')
    print(f"\nComparison figures saved as: {output_filename}.pdf")
    cursor = mplcursors.cursor([ax1, ax2, ax3, ax4, ax5, ax6], hover=True)

    plt.show(block=False)
    plt.pause(500)
    #plt.close()

if __name__ == "__main__":
    main()