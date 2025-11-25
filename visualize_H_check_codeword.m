%% Author: Njini Nathan Fofeyin %%

B = load("base_matrix.txt");
H = obtain_H(B);
spy(H, 15);
xlabel("$n$-Index", Interpreter="latex", FontSize=12)
ylabel("$m$-Index", Interpreter="latex", FontSize=12)
title("$1152 \times 2304$ Parity Check Matrix $\mathbf{H}$", Interpreter="latex", FontSize=14)