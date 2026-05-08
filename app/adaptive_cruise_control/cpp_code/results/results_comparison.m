% Example main script to load and plot data
    bipm_filename =  "../final/260216_1119_G2O_ISPD_gn_GN_12_3.bin"
    al_filename = "260217_1315_AMPL_0_3.bin"  

bipm_data = load_bin_data(bipm_filename);
al_data = load_bin_data(al_filename);

figure;
subplot(1,3,1);
plot(bipm_data.v_p_signal, 'b-', 'DisplayName', 'BIPM v_p'); hold on;
plot(bipm_data.v_h_signal, 'b--', 'DisplayName', 'BIPM v_h');
plot(al_data.v_h_signal, 'r-', 'DisplayName', 'AL v_h');
ylabel('Velocity (m/s)'); legend; grid on;

subplot(1,3,2);
plot(bipm_data.num_iterations_signal, 'b-', 'DisplayName', 'BIPM Iterations'); hold on;
plot(al_data.num_iterations_signal, 'r-', 'DisplayName', 'AL Iterations');
ylabel('Iterations'); legend; grid on;

subplot(1,3,3);
plot(bipm_data.d_h_signal, 'b-', 'DisplayName', 'BIPM Distance'); hold on;
plot(al_data.d_h_signal, 'r-', 'DisplayName', 'AL Distance');
ylabel('Inter-distance (m)'); legend; grid on;

sgtitle('MPC Performance Comparison: BIPM vs AL');




