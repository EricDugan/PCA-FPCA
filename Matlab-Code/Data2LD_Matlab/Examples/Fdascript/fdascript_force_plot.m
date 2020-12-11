function fdascript_force_plot(Aknots, Acoef, Aname, vertrng)
nAknots = length(Aknots);
k = 1;
phdl = plot([Aknots(k),Aknots(k+1)], ...
            [Acoef(k) ,Acoef(k)], 'b-');
set(phdl, 'LineWidth', 2)
axis([0,230,vertrng])
hold on
phdl = plot([Aknots(k+1),Aknots(k+1)], ...
            [Acoef(k),   Acoef(k+1)], 'b-');
set(phdl, 'LineWidth', 2)
for k=2:nAknots-2
    phdl = plot([Aknots(k),  Aknots(k+1)], ...
                [Acoef(k),   Acoef(k)],   'b-');
    set(phdl, 'LineWidth', 2)
    phdl = plot([Aknots(k+1),Aknots(k+1)], ...
                [Acoef(k),   Acoef(k+1)], 'b-');
    set(phdl, 'LineWidth', 2)
end
phdl = plot([Aknots(nAknots-1),Aknots(nAknots)], ...
            [Acoef(nAknots-1), Acoef(nAknots-1)], 'b-', ...
            [Aknots(1),Aknots(nAknots)], [0,0], 'b--');
set(phdl, 'LineWidth', 2)
hold off
ylabel(['\fontsize{16} ',Aname])
