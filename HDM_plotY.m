function [] = HDM_plotY(y)

p = HDM_getParameters();

titles = {'n_{excitation}', 'n_{inhibition}', 'vaso', 'f_{arteriole}', 'f', '-', 'v', '-', 'q', '-', 'signal'};
for f = 1:length(titles)
    if strcmpi(titles{f},'-')
        continue;
    end
    figure;
    for d = 1:p.D
        subplot(p.D,1,d);
        if ~any(strcmpi(titles{f}, {'f','v','q'}))
            plot(y(d+(f-1)*p.D,:));
        else
            plot(y(d+(f-1)*p.D,:), 'color','magenta');
            plot(y(d+(f-1)*p.D+p.D,:), 'color','blue');
        end
        title(titles{f});
    end
end

end

