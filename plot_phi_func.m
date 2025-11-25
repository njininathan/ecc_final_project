%% Author: Njini Nathan Fofeyin %%
x = 0:0.01:6;
y = log ((exp(x) + 1)./(exp(x) - 1));

plot( x, y, LineWidth=2);
xlabel("$ x $", Interpreter="latex",FontSize=14)
ylabel("$\phi (x)$", Interpreter="latex", FontSize=14, Rotation=0)
hold on;
plot(x,x, LineStyle="--", Color="#800080");
gtext("45°")
hold off;
grid;