function aux_func_plot_circle(x, y, r, id_figure, opt_plot)
figure(id_figure);
hold on
if r > 0
    ds = 0:0.01:2*pi;
    x_edge = x + r * cos(ds);
    y_edge = y + r * sin(ds);
    plot(x_edge, y_edge, opt_plot);
end
% plot(x,y,'rd');
% if r > 0
%     xlim([min(x_edge) max(x_edge)]);
%     ylim([min(y_edge) max(y_edge)]);
% end
end