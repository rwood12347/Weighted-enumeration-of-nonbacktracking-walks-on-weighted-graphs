%Code for the normalised NBT/Normalised Katz centrality graphs at t
%=0.95/rho for binarised and non-binarised graphs.
clear all


load('D:\Fauci paper redo\Static Networks\Hypergraph-projection-cc\Matrices\M.mat');
M_struc = load('D:\Fauci paper redo\Static Networks\Hypergraph-projection-cc\Matrices\binarised M.mat');
bin_M = M_struc.M;
load('D:\Fauci paper redo\Static Networks\Hypergraph-projection-cc\Matrices\TargetR.mat');
load('D:\Fauci paper redo\Static Networks\Hypergraph-projection-cc\Matrices\SourceL.mat');
load('D:\Fauci paper redo\Static Networks\Hypergraph-projection-cc\Matrices\A.mat');
binarised_a = load('D:\Fauci paper redo\Static Networks\Hypergraph-projection-cc\Matrices\binarised A.mat');
bin_A = binarised_a.A;
load('D:\Fauci paper redo\Static Networks\Hypergraph-projection-cc\Matrices\sqrtS.mat');
load('D:\Fauci paper redo\Static Networks\Hypergraph-projection-cc\Matrices\rhos_KatzNBT.mat');
binarised_rhos = load('D:\Fauci paper redo\Static Networks\Hypergraph-projection-cc\Matrices\binarised rhos_KatzNBT.mat');
bin_rho_katz = binarised_rhos.both_rhos(1);
bin_rho_nbt = binarised_rhos.both_rhos(2);
load('D:\Fauci paper redo\Static Networks\Hypergraph-projection-cc\Matrices\labels.mat');
tic
display('Files loaded')
n = size(A,1);
m = size(M,1);


%Make the ranges of attenuation factors

Katz_t = 0.95/abs(both_rhos(1));
Bin_Katz_t = 0.95/bin_rho_katz;

NBT_t =  0.95/abs(both_rhos(2));
Bin_NBT_t = 0.95/bin_rho_nbt;

display("t's loaded.")

%The following code computes the centrality measure for the above range.
%


resolvent = inv(speye(size(M)) - NBT_t*M);
disp('Resolvent made.')
x_NBT = (speye(n) + NBT_t*SourceL.'*sqrtS*resolvent*sqrtS*TargetR)*ones(n,1);
disp('NBT Centrality Computed.')
clear resolvent
resolvent = inv(speye(size(M)) - Bin_NBT_t*bin_M);
bin_x_NBT = (speye(n) + Bin_NBT_t*SourceL.'*resolvent*TargetR)*ones(n,1);
%The following loop calculates Katz centrality with attentuation factor t
%for the time evolving network.
x_katz = Katz_method(A, Katz_t, ones(n,1) );
bin_x_katz = Katz_method(bin_A, Bin_Katz_t, ones(n,1) );
%normalise the centrality vectors

x_katz = x_katz / norm(x_katz,2);
x_NBT = x_NBT/norm(x_NBT,2);

    
bin_x_katz = bin_x_katz / norm(bin_x_katz,2);
bin_x_NBT = bin_x_NBT/norm(bin_x_NBT,2);


% 
% % NBT Graph
% 
% figure 
% 
% 
% NBTplot = plot(graph(A),  'Layout',"force", 'NodeCData', x_NBT);
% NBTplot.MarkerSize = x_NBT;
% colorbar
% 
% %Label only the most prominent nodes.
% for i=1:size(A,2)
% if  x_NBT(i)  < 0.1*max(x_NBT)
%     newlabels{i} = " ";
% else
%     newlabels{i} = labels(i);
% end
% end
% 
% labelnode(NBTplot,1:size(A,2), string(newlabels))
% 
% 
% 
% % Katz Graph
% 
% figure 
% 
% 
%  Katzplot = plot(graph(A),  'Layout',"force", 'NodeCData', x_katz)
%  Katzplot.MarkerSize = x_katz
% colorbar
% 
% for i=1:size(A,2)
% if  x_katz(i)  < 0.1*max(x_katz)
%     newlabels{i} = " ";
% else
%     newlabels{i} = labels(i);
% end
% end
% 
% labelnode(Katzplot,1:size(A,2), string(newlabels))
% 
% 
% 
% 
% 
no_top = 5

[~,Katz_top_indices] = sort(x_katz, 'descend');
Katz_top_indices = Katz_top_indices(1:no_top);


[~,NBT_top_indices] = sort(x_NBT, 'descend');
NBT_top_indices = NBT_top_indices(1:no_top);




union_top_indices = Katz_top_indices;
weighted_cut_off = min(x_katz(union_top_indices))



[~,bin_Katz_top_indices] = sort(bin_x_katz, 'descend');
bin_Katz_top_indices = bin_Katz_top_indices(1:no_top);


[~,bin_NBT_top_indices] = sort(bin_x_NBT, 'descend');
bin_NBT_top_indices = bin_NBT_top_indices(1:no_top);

bin_union_top_indices = union(bin_NBT_top_indices,bin_Katz_top_indices);
bin_cut_off = min(min(bin_x_NBT(bin_union_top_indices)), min(bin_x_katz(bin_union_top_indices)))

for i = 1:size(labels,2)
    if (max(x_NBT(i),x_katz(i)) < weighted_cut_off)
weighted_labels{i} = " ";
    else
        weighted_labels{i} = labels(i);
    end

end
for i = 1:size(labels,2)
    if (min(bin_x_NBT(i),bin_x_katz(i)) < bin_cut_off)
bin_labels{i} = " ";
    else
        bin_labels{i} = labels(i);
    end

end

% 
% I = eye(20);
% % P = [I(:,1), I(:,2), I(:,11), I(:,3), I(:,4), I(:,12), I(:,5), I(:,13), I(:,14), I(:,15), I(:,6), I(:,7), I(:,8), I(:,9), I(:,10), I(:,16), I(:,17), I(:,18), I(:,19), I(:,20)];
% % union_top_indices = P*union_top_indices;
% 
% union_top_values = [x_katz(union_top_indices), x_NBT(union_top_indices)];
% union_labels = labels(union_top_indices);
% 


figure
subplot(1,2,2)
scatterplot = scatter(x_katz, x_NBT, 'filled')
text(x_katz+0.005*ones(891,1), x_NBT+0.005*ones(891,1), weighted_labels,'FontSize', 10)
ax = gca
ax.YAxis.Label.String = 'Normalized Non-backtracking Katz Centrality Value';
ax.XAxis.Label.String = 'Normalized Katz Centrality Value';
ax.Title.String = 'Weighted scatter';
ax.XLim = [0,0.55]
ax.YLim = [0,0.55]
hold on
plot(0:0.05:0.55, 0:0.05:0.55, 'LineStyle', '-.')
hold off

subplot(1,2,1)
bin_scatterplot = scatter(bin_x_katz, bin_x_NBT, 'filled')
text(bin_x_katz+0.003*ones(891,1), bin_x_NBT+0.001*ones(891,1), bin_labels,'FontSize', 10)
ax = gca
ax.YAxis.Label.String = 'Normalized Non-backtracking Katz Centrality Value';
ax.XAxis.Label.String = 'Normalized Katz Centrality Value';
ax.Title.String = 'Binarised scatter';
ax.XLim = [0, 0.35]
ax.YLim = [0, 0.35]
hold on
plot(0:0.05:0.35, 0:0.05:0.35, 'LineStyle', '-.')
hold off

toc