function [] = HDM_plotD(tit, mat1, mat2, legend1, legend2 )

p = HDM_getParameters();

figure;

for d = 1:p.D
    subplot(p.D,1,d);
    plot(mat1(d,:));
    hold on;
    if exist('mat2', 'var')
        plot(mat2(d,:));
        if exist('legend2', 'var') && d==1
            title(tit);
            legend(legend1, legend2);
        end
    end
end


end

