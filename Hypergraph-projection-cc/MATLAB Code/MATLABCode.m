clear all

load('D:\Fauci paper redo\Static Networks\Hypergraph-projection-cc\Matrices\M.mat');
load('D:\Fauci paper redo\Static Networks\Hypergraph-projection-cc\Matrices\TargetR.mat');
load('D:\Fauci paper redo\Static Networks\Hypergraph-projection-cc\Matrices\SourceL.mat');
load('D:\Fauci paper redo\Static Networks\Hypergraph-projection-cc\Matrices\A.mat');
load('D:\Fauci paper redo\Static Networks\Hypergraph-projection-cc\Matrices\sqrtS.mat');
load('D:\Fauci paper redo\Static Networks\Hypergraph-projection-cc\Matrices\rhos_KatzNBT.mat');
load('D:\Fauci paper redo\Static Networks\Hypergraph-projection-cc\Matrices\labels.mat');
load('Katzcent.mat');
load('NBTcent.mat');
n = size(A,1);
m = size(M,1);

%Make the ranges of attenuation factors

Katz_t_range = [0, 0.5/abs(both_rhos(1)), 0.75/abs(both_rhos(1)), 0.95/abs(both_rhos(1))];
NBT_t_range =[0, 0.5/abs(both_rhos(2)), 0.75/abs(both_rhos(2)), 0.95/abs(both_rhos(2))];


no_steps = size(NBT_t_range,2);

%The following code computes the centrality measure for the above range.
%

%loop over the above ranges
% for k = 1:no_steps
%     k
%     Katz_t = Katz_t_range(k);
%     NBT_t = NBT_t_range(k);
%     resolvent = inv(speye(size(M)) - NBT_t*M);
%     disp('resolvent made')
%     x_NBT = (eye(n) + NBT_t*SourceL.'*sqrtS*resolvent*sqrtS*TargetR)*ones(n,1);
%     disp('cent done ')
%     
% %The following loop calculates Katz centrality with attentuation factor t
% %for the time evolving network.
% 
% 
% x_katz = Katz_method(A, Katz_t, ones(n,1) );
% 
% %normalise the centrality vector
% 
% x_katz = x_katz / norm(x_katz,2);
% x_NBT = x_NBT/norm(x_NBT,2);
% 
%     
% Katz(:,k) = x_katz;
% NBT(:,k) = x_NBT;
% end


k = 2





[~,Katz_top_indices] = sort(Katz(:,k), 'descend');
Katz_top_indices = Katz_top_indices(1:10);


[~,NBT_top_indices] = sort(NBT(:,k), 'descend');
NBT_top_indices = NBT_top_indices(1:10);

union_top_indices = union(NBT_top_indices,Katz_top_indices);


I = eye(20);
% P = [I(:,1), I(:,2), I(:,11), I(:,3), I(:,4), I(:,12), I(:,5), I(:,13), I(:,14), I(:,15), I(:,6), I(:,7), I(:,8), I(:,9), I(:,10), I(:,16), I(:,17), I(:,18), I(:,19), I(:,20)];
% union_top_indices = P*union_top_indices;

union_top_values = [NBT(union_top_indices,k), Katz(union_top_indices,k)];
union_labels = labels(union_top_indices);



subplot(1,2,1)
barchart = bar(union_top_values)


barchart = bar(categorical(union_labels), union_top_values, 'LineStyle', "--", 'LineWidth', 1.5, 'BarWidth', 1)
barchart(1).LineStyle = ":";
legend('NBT Katz','Katz')
ax = gca;
ax.FontSize = 24;
ax.Legend.FontSize = 24;
ax.Legend.LineWidth = 1.5;
ax.Title.String = sprintf("Katz centrality values with attenuation factor t = 0.5/ρ");
ax.Title.FontSize = 24;
ax.YAxis.Label.String = 'Normalized centrality value';
ax.Title.Position(2) = ax.Title.Position(2)*1.02
%ax.XAxis.Label.String = 'Union of 10 most highly ranked nodes for Katz and NBT Katz';


k = 4

[~,Katz_top_indices] = sort(Katz(:,k), 'descend');
Katz_top_indices = Katz_top_indices(1:10);


[~,NBT_top_indices] = sort(NBT(:,k), 'descend');
NBT_top_indices = NBT_top_indices(1:10);

union_top_indices = union(NBT_top_indices,Katz_top_indices);


I = eye(20);
% P = [I(:,1), I(:,2), I(:,11), I(:,3), I(:,4), I(:,12), I(:,5), I(:,13), I(:,14), I(:,15), I(:,6), I(:,7), I(:,8), I(:,9), I(:,10), I(:,16), I(:,17), I(:,18), I(:,19), I(:,20)];
% union_top_indices = P*union_top_indices;

union_top_values = [NBT(union_top_indices,k), Katz(union_top_indices,k)];
union_labels = labels(union_top_indices);




subplot(1,2,2)
barchart = bar(union_top_values)


barchart = bar(categorical(union_labels), union_top_values, 'LineStyle', "--", 'LineWidth', 1.5, 'BarWidth', 1)
barchart(1).LineStyle = ":";
legend('NBT Katz','Katz')
ax = gca;
ax.FontSize = 24;
ax.Legend.FontSize = 24;
ax.Legend.LineWidth = 1.5;
ax.Title.String = sprintf("Katz centrality values with attenuation factor t = 0.95/ρ");
ax.Title.FontSize = 24;
ax.Title.Position(2) = ax.Title.Position(2)*1.02
ax.YAxis.Label.String = 'Normalized centrality value';
%ax.XAxis.Label.String = 'Union of 10 most highly ranked nodes for Katz and NBT Katz';
