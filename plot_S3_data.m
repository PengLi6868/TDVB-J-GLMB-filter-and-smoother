clear all
clc,close all

A{1} = [11.55	11.16	8.903	9.033
    13.3	12.27	10.97	9.784
    14.98	13.31	12.38	10.51
    17.68	14.4	14.48	11.67
    19.59	15.68	16.7	13.35];
A{2} = [11.55	11.16	8.903	9.033
    14.52	13.65	11.82	10.97
    15.35	14.08	12.73	12.32
    18.75	16.87	15.55	15.48
    22.98	18.89	21.55	17.63];
A{3} = [19.88	16.36	16.26	12.27
    23.65	17.93	19.41	12.78
    26.59	19.1	21.77	14.6
    33.27	21.78	29.4	17.13
    33.05	23.26	29.08	19.3];
A{4} = [19.88	16.36	16.26	12.27
    20.83	17.53	15.84	13.72
    24.67	20.31	21.01	18.23
    27.28	22.8	25.18	20.06
    30.9	25.83	27.88	24.33];

N = length(A);

c = {'b-*','r-o','g->','k-s'};
for i = 1:N
    figure(10+i)
    hold on
    box on
    for j = 1:4
        plot([1:5],A{i}(:,j)',c{j})
    end
    legend('J-GLMB','TDVB-J-GLMB','J-GLMB with smoothing','TDVB-J-GLMB with smoothing')
    if i == 1 || i == 3
        title('Performance under different cluter rate (p_D = 0.95)')
        xlabel('Clutter rate'),ylabel('Average OSPA')
        set(gca,'XTick',1:5);set(gca,'XTicklabel',{'20','40','60','80','100'})
    else
        title('Performance under different detection probability (¦Ã = 20)')
        xlabel('Detection probability'),ylabel('Average OSPA')
        set(gca,'XTick',1:5);set(gca,'XTicklabel',{'0.95','0.90','0.85','0.80','0.75'})
    end
end
