%% Author: Njini Nathan Fofeyin %%
Eb_N0_dB = [0 0.2 0.4 0.6 0.8 1 1.2 1.4 1.6 1.8 2 2.2];

log_SPA_BER = [1.230e-001 1.109e-001 9.672e-002 7.907e-002 5.682e-002 3.209e-002 1.264e-002 3.440e-003 6.371e-004 6.226e-005 5.660e-006 6.097e-007];
minsum_BER = [1.752e-01 1.658e-01 1.551e-01 1.435e-01 1.304e-01 1.136e-01 8.966e-02 5.737e-02 2.602e-02 6.601e-003  1.077e-003 8.950e-005];

figure;
semilogy(Eb_N0_dB, log_SPA_BER, '-o', LineWidth=2, MarkerSize=6)
hold on;
semilogy(Eb_N0_dB, minsum_BER, '-o', LineWidth=2, MarkerSize=6, Color="#800080")
grid;
xlabel("EbN0 (dB)", Interpreter="latex")
ylabel("BER", Interpreter="latex")
title("Log SPA vs minsum Comparison", Interpreter="latex", FontSize=14)
legend("Log SPA", "minsum", Interpreter="latex", FontSize= 12)
hold off;

figure;
semilogy(Eb_N0_dB, log_SPA_BER, '-o', LineWidth=2, MarkerSize=6)
legend("LDPC rate = 0.5", Interpreter="latex", FontSize= 12)
title("LDPC Log SPA Simulation BER vs $\frac{E_b}{N_0}$", Interpreter="latex", FontSize=14)
grid;