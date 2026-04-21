import time  
import os
import sys
import numpy as np
import argparse
import matplotlib.pyplot as plt
from data_saver.data_loader import DataLoader
from datetime import datetime

def safe_get(results, key, default="N/A"):
    """Safely get a value from results dictionary with default fallback"""
    return results.get(key, default)

def safe_get_float(results, key, default=0.0):
    """Safely get a float value with default fallback"""
    try:
        return float(results.get(key, default))
    except (ValueError, TypeError):
        return default

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Load and plot MPC data.")
    
    # --data flag (defaulting to your specific test file)
    parser.add_argument('--data', type=str, 
                        default="260221_0747_G2O_ISPD_gn_GN_13_3.bin",
                        help="The filename of the data to load")
    
    # --dir flag (defaulting to the 'data' folder relative to the script)
    parser.add_argument('--dir', type=str, 
                        default=os.path.join(os.path.dirname(__file__), 'data'),
                        help="The directory where the data file is located")

    args = parser.parse_args()

    # Use the arguments
    filename = args.data
    data_dir = args.dir
    
    # Construct the full path
    data_file = os.path.join(data_dir, filename)
    
    print(f"Loading data from: {data_file}")
    
    # Check if file exists to prevent DataLoader crash
    if not os.path.exists(data_file):
        print(f"Error: File not found at {data_file}")
        sys.exit(1)
    
    print_results = False
    # Load all data first
    with DataLoader(data_file) as loader:
        results = loader.load_all()
        if print_results:
            print(f"Loaded {len(results)} variables from {filename}")
            
            # Print available variables
            print("=" * 50)
            print("Available variables:")
            print("-" * 50)
            for i, name in enumerate(results.keys(), 1):
                print(f"{i}. {name}")
            print("=" * 50)

    # ======================================================================
    # Safely load all parameters with defaults
    # ======================================================================
    file_name = safe_get(results, 'File_name')
    solverType = safe_get(results, 'Alg_unconstrainted')
    AlgEqType = safe_get(results, 'Alg_equality')
    print(f"-------------------------------------------AlgEqType: {AlgEqType}")
    AlgIneqType = safe_get(results, 'Alg_inequality')
    terminationCriterionStr = safe_get(results, 'Termination_criteria:')
    linearSolverType = int(safe_get(results, 'linearSolverType', -1))
    mpcHorizon = int(safe_get(results, 'mpcHorizon', -1))
    
    # Signals (initialize as empty arrays if missing)
    v_h_signal = np.array(safe_get(results, 'v_h_signal', []))
    d_h_signal = np.array(safe_get(results, 'd_h_signal', []))
    v_p_signal = np.array(safe_get(results, 'v_p_signal', []))
    a_p_signal = np.array(safe_get(results, 'a_p_signal', []))
    num_iterations_signal = np.array(safe_get(results, 'num_iterations_signal', []))
    cal_time_signal = np.array(safe_get(results, 'cal_time_signal', []))
    u_input_signal = np.array(safe_get(results, 'u_input_signal', []))
      


    # print("u_input_signal:", u_input_signal)
    # print("cal_time_signal:", cal_time_signal)
    # print("v_h_signal:", v_h_signal)
    # print("num_iterations_signal:", num_iterations_signal)
    # AL-specific parameters
    settingEq0 = safe_get_float(results, 'init_lagrange_multiplier_eq',-1)
    settingEq1 = safe_get_float(results, 'init_rho_eq',-1)
    settingEq2 = safe_get_float(results, 'rho_min_eq',-1)
    settingEq3 = safe_get_float(results, 'rho_max_eq',-1)
    settingEq4 = safe_get_float(results, 'rho_update_factor_eq',-1)
    settingEq5 = safe_get_float(results, 'max_num_inner_iterations',-1)
    
    settingIneq0 = safe_get_float(results, 'init_lagrange_multiplier_ineq',-1)
    settingIneq1 = safe_get_float(results, 'init_rho_ineq',-1)
    settingIneq2 = safe_get_float(results, 'rho_min_ineq',-1)
    settingIneq3 = safe_get_float(results, 'rho_max_ineq',-1)
    settingIneq4 = safe_get_float(results, 'rho_upate_factor_ineq',-1)
    
    # BIPM-specific parameters
    t_final = safe_get_float(results, 't_final',-1)
    t_init_bipm = safe_get_float(results, 't_init',-1)
    nu_update = safe_get_float(results, 'nu_update',-1)
    
    #ISPD-specific parameters
    linear_system_strategy = safe_get_float(results,'linear_system_strategy', -1)
    step_size_strategy = safe_get_float(results,'step_size_strategy', -1)
    ineq_backtracking_step_min = safe_get_float(results,'ineq_backtracking_step_min', -1)
    aux_backtracking_step_min = safe_get_float(results,'aux_backtracking_step_min', -1)
    
    t_0_strategy = safe_get_float(results,'_init_kappa_strategy', -1)
    t_init_ispd = safe_get_float(results,'t_init', -1)
 

    update_kappa_strategy = safe_get_float(results,'update_kappa_strategy', -1)
    nu = safe_get_float(results,'nu_update', -1)
    t_final = safe_get_float(results,'t_final', -1)
    gamma = safe_get_float(results,'gamma', -1)
    limit_kappa_final = safe_get_float(results,'limit_kappa_final', -1);
    
    
    
    init_slack_strategy = safe_get_float(results,'init_slack_strategy', -1)
    s_0 = safe_get_float(results,'slack_variable_initial_ineq', -1)
    
    
    init_lagrange_strategy = safe_get_float(results,'init_lagrange_strategy', -1)
    lambda_0 = safe_get_float(results,'lagrange_multiplier_initial_ineq', -1)
    
    error_pred_strategy = safe_get_float(results,'ineq_prediction_strategy', -1);
    ineq_prediction_step_strategy = safe_get_float(results,'ineq_prediction_step_strategy', -1)
    error_pred_param = safe_get_float(results,'ineq_prediction_param', -1);
    
    keep_aux_positive_strategy = safe_get_float(results,'keep_aux_positive_strategy', -1)
    aux_scaling_factor = safe_get_float(results,'aux_scaling_factor', -1)
    aux_correction_value = safe_get_float(results,'aux_correction_value', -1)
    
    
    feas_tol = safe_get_float(results,'feas_tol', -1); 
     

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
    if len(v_h_signal) > 0 or len(v_p_signal) > 0:
        plt.figure(figsize=(15, 12))
        plt.suptitle('MPC Performance Analysis', y=0.98, fontsize=16)

        # plot index form zero to first zero in calcuation time sginal 
        if len(cal_time_signal) > 0:
            cal_time_signal = cal_time_signal[cal_time_signal > 0]
            num_iterations_signal = num_iterations_signal[:len(cal_time_signal)]
            v_h_signal = v_h_signal[:len(cal_time_signal)]
            v_p_signal = v_p_signal[:len(cal_time_signal)]
            d_h_signal = d_h_signal[:len(cal_time_signal)]
            u_input_signal = u_input_signal[:len(cal_time_signal)]
        # Subplot 1: v_p and v_h
        ax1 = plt.subplot(3, 2, 1)
        if len(v_p_signal) > 0:
            ax1.plot(v_p_signal, 'b-', label='v_p (preceding vehicle)')
        if len(v_h_signal) > 0:
            ax1.plot(v_h_signal, 'r-', label='v_h (host vehicle)')
        ax1.set_ylabel('Velocity (m/s)')
        ax1.legend()
        ax1.grid(True)

        # Subplot 2: Iterations with mean/std
        ax2 = plt.subplot(3, 2, 2)
        if len(num_iterations_signal) > 0:
            iter_mean = np.nanmean(num_iterations_signal)
            iter_std = np.nanstd(num_iterations_signal)
            ax2.plot(num_iterations_signal, 'y-', label='Iterations')
            ax2.plot([iter_mean]*len(num_iterations_signal), 'r-', label='Mean')
            ax2.plot([iter_mean + iter_std]*len(num_iterations_signal), 'r--')
            ax2.plot([iter_mean - iter_std]*len(num_iterations_signal), 'r--')
        ax2.set_ylabel('Iterations')
        ax2.legend()
        ax2.grid(True)

        # Subplot 3: Distance
        ax3 = plt.subplot(3, 2, 3)
        if len(d_h_signal) > 0:
            ax3.plot(d_h_signal, 'g-')
        ax3.set_ylabel('Inter-distance (m)')
        ax3.grid(True)

        # Subplot 4: Computation time with mean/std
        ax4 = plt.subplot(3, 2, 4)
        if len(cal_time_signal) > 0:
            time_mean = np.nanmean(cal_time_signal)
            time_std = np.nanstd(cal_time_signal)
            ax4.plot(cal_time_signal, 'c-', label='Comp Time')
            ax4.plot([time_mean]*len(cal_time_signal), 'r-', label='Mean')
            ax4.plot([time_mean + time_std]*len(cal_time_signal), 'r--')
            ax4.plot([time_mean - time_std]*len(cal_time_signal), 'r--')
        ax4.set_ylabel('Time (ms)')
        ax4.legend()
        ax4.grid(True)

        # Subplot 5: Control input
        ax5 = plt.subplot(3, 2, 5)
        if len(u_input_signal) > 0:
            ax5.plot(u_input_signal, 'm-')
        ax5.set_xlabel('Time steps')
        ax5.set_ylabel('Control Input (N)')
        ax5.grid(True)

        # Subplot 6: Time per iteration with mean/std
        ax6 = plt.subplot(3, 2, 6)
        if len(cal_time_signal) > 0 and len(num_iterations_signal) > 0:
            cal_time_per_iter = np.array(cal_time_signal) / np.array(num_iterations_signal)
            ct_mean = np.nanmean(cal_time_per_iter)
            ct_std = np.nanstd(cal_time_per_iter)
            ax6.plot(cal_time_per_iter, 'g-', label='Time/Iter')
            ax6.plot([ct_mean]*len(cal_time_per_iter), 'r-', label='Mean')
            ax6.plot([ct_mean + ct_std]*len(cal_time_per_iter), 'r--')
            ax6.plot([ct_mean - ct_std]*len(cal_time_per_iter), 'r--')
            if print_results:
                print(f"Mean time per iteration: {ct_mean:.2f} ms")
        ax6.set_xlabel('Time steps')
        ax6.set_ylabel('Time/Iter (ms)')
        ax6.legend()
        ax6.grid(True)

        # ======================================================================
        # Text Configuration Section
        # ======================================================================
        # Algorithm configuration text
        text_lines = [
        f"Algorithm Configuration:",
        f"- Inequality: {AlgIneqType} | Equality: {AlgEqType} | Unconstrained: {solverType}  | Linear Solver: {solver_types.get(linearSolverType, 'Unknown')} (Type {linearSolverType}) | Horizon: {mpcHorizon} | Termination: {terminationCriterionStr}"
        ]

        # Add AL parameters if available
        if AlgIneqType == 'AL':
            text_lines.extend([
                f"- AL Ineq Params: λ_init: {settingIneq0:.0e} | ρ_init: {settingIneq1} | ρ_min: {settingIneq2} | ρ_max: {settingIneq3:.0e} | ρ_upate_factor: {settingIneq4}"            
            ])
             
        # Add BIPM parameters if available and relevant
        if AlgIneqType == 'BIPM':
            text_lines.extend([
                f"- BIPM Params: t_final: {t_final:.1f} | t_0: {t_init_bipm:.1f} | ν_update: {nu_update:.1f}"                   
             ])
            
        if AlgIneqType == 'ISPD':   
            text_lines.extend([
                f"- ISPD Params: lin_sys_s: {linear_system_strategy:.0f} | step_size_s: {step_size_strategy:.0f} | ineq_step_min: {ineq_backtracking_step_min:.0e} | aux_step_min: {aux_backtracking_step_min:.1f} | t_0_s: {t_0_strategy:.0f} | t_0: {t_init_ispd:.0e} | t_update_s: {update_kappa_strategy:.0f} | nu: {nu:.1f} | t_final: {t_final:.0f}"
        ])  
            text_lines.extend([
                f"- Cont.: | gamma: {gamma:.2f} | limit_kappa_final: {limit_kappa_final:.1f} | slack_init_s: {init_slack_strategy:.0f} | s_0: {s_0:.0e} | lambda_init_s: {init_lagrange_strategy:.0f} | lambda_0: {lambda_0:.0e} | error_pred_s: {error_pred_strategy:.0f} | error_pred_step_s: {ineq_prediction_step_strategy:.0f}"
        ])  
            text_lines.extend([
                f"- Cont.: | error_pred_param: {error_pred_param:.1f} | keep_aux_pos_s: {keep_aux_positive_strategy:.0f} | aux_scaling: {aux_scaling_factor:.1f} | aux_corr: {aux_correction_value:.1f} | error_pred_param: {error_pred_param:.1f} | feas_tol: {feas_tol:.3f}"
        ])  
            
            
            
            

        # Add AL inequality parameters if available and relevant
        try:
            if results['Alg_equality'] == 'al':
                text_lines.extend([
                f"- AL Eq Params: λ_init: {settingEq0:.1e} | ρ_init: {settingEq1:.1f} | ρ_min: {settingEq2:.1f} | ρ_max: {settingEq3:.1e} | ρ_upate_factor: {settingEq4:.1f} | max_num_inner_iterations: {settingEq5:.1f}"
            ])
            else:   
                text_lines.extend([
                    f"- Eq Params: max_num_inner_iterations: {settingEq5:.1f}"
                ])
        except KeyError:
            text_lines.extend([
                f"- Eq Params: max_num_inner_iterations: {settingEq5:.1f}"
            ])
       

        # Add statistics if available
        if len(num_iterations_signal) > 0:
            text_lines.append(
                f"- Iterations: mean={np.nanmean(num_iterations_signal):.1f} | min={np.nanmin(num_iterations_signal)} | max={np.nanmax(num_iterations_signal)} | σ={np.nanstd(num_iterations_signal):.1f})"
            )

        if len(cal_time_signal) > 0:
            text_lines.append(
                f"- Comp Time: mean={np.nanmean(cal_time_signal):.1f}ms | min={np.nanmin(cal_time_signal):.1f} | max={np.nanmax(cal_time_signal):.1f} | σ={np.nanstd(cal_time_signal):.1f})"
            )

        if len(cal_time_signal) > 0 and len(num_iterations_signal) > 0:
            text_lines.append(f"- Time/Iter: mean={ct_mean:.2f}ms | min={np.nanmin(cal_time_per_iter):.1f} | max={np.nanmax(cal_time_per_iter):.1f} || (σ={ct_std:.2f})"
            )

        plt.figtext(0.5, 0.02, '\n'.join(text_lines), 
                   ha='center', fontsize=9, 
                   bbox=dict(facecolor='lightgray', alpha=0.5))

        plt.tight_layout()
        plt.subplots_adjust(bottom=0.2)

        # Save figure
        #timestamp = datetime.now().strftime("%y%m%d_%H%M")
        #output_filename = f"{timestamp}_{AlgIneqType}_{AlgEqType}_{solverType}_{linearSolverType}_{mpcHorizon}"
        base_name = os.path.splitext(filename)[0]
        output_path = os.path.join(args.dir, base_name)
        save_format = 1  # Replace with the appropriate logic to determine the format
        if save_format == 1:
            full_output = f"{output_path}.pdf"
            plt.savefig(full_output, format='pdf', bbox_inches='tight')
            print(f"\nFigure saved in data directory: {full_output}")
        elif save_format == 2:
            full_output = f"{output_path}.png"
            plt.savefig(full_output, dpi=300, bbox_inches='tight')
            print(f"\nFigure saved in data directory: {full_output}")
        
        plt.show(block=False)  # Non-blocking
        plt.pause(15)  # Show for 3 seconds
        plt.close()
    
    else:
        print("the size of v_h_signal is:", len(v_h_signal), "and the size of v_p_signal is:", len(v_p_signal))
        
        print("\nNo signal data available for plotting")


if __name__ == "__main__":
    main()
 
