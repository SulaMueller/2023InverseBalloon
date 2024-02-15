function [] = HDM_plotD(D, tit, mat1, mat2, legend1, legend2)
%% function [] = HDM_plotD(D, tit, mat1, mat2, legend1, legend2)
% funtion to plot 1 or 2 functions with D (number of depth levels) subplots

figure;

for d = 1:D
    subplot(D,1,d);
    plot(mat1(d,:),'LineWidth',4);
    hold on;
    if exist('mat2', 'var')
        plot(mat2(d,:),'LineWidth',4);
        if exist('legend2', 'var') && d==1
            legend(legend1, legend2);
        end
    end
    if d==1 
        title(tit);
    end
end

end

