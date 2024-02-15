function [] = HDM_plotY(D, y)
%% function [] = HDM_plotY(D, y)
% function to plot all info contained in Y (state matrix for all state
% variables at all time points; for definition view HDM_solveForward)

%% define structure of Y
% '-', if preceding variable needs two slots (venule/vein)
titles = {'n_{excitation}', 'n_{inhibition}', 'vaso', 'f_{arteriole}', 'f', '-', 'v', '-', 'q', '-', 'signal'};

%% plot
for f = 1:length(titles)
    if strcmpi(titles{f},'-')  % has been plotted in preceding loop
        continue;
    end
    figure;
    for d = 1:D
        subplot(D,1,d);
        if ~any(strcmpi(titles{f}, {'f','v','q'}))
            plot(y(d+(f-1)*D,:),'LineWidth',4);
        else
            plot(y(d+(f-1)*D,:), 'color','magenta','LineWidth',4);  % venule
            plot(y(d+(f-1)*D+D,:), 'color','blue','LineWidth',4);  % vein
        end
        if d==1
            title(titles{f});
        end
    end
end

end

