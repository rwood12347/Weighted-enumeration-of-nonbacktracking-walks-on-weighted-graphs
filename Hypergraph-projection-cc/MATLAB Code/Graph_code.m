load('D:\Fauci paper redo\Static Networks\Hypergraph-projection-cc\Matrices\A.mat');
load('Katzcent.mat');
load('NBTcent.mat');
load('D:\Fauci paper redo\Static Networks\Hypergraph-projection-cc\Matrices\labels.mat');
%NBT Graph

figure 


NBTplot = plot(graph(A),  'Layout',"force", 'NodeCData', NBT(:,k));
NBTplot.MarkerSize =   ( NBT(:,k)/norm( NBT(:,k),2));
colorbar

for i=1:size(A,2)
if  NBT(i,k)  < 0.1*max( ( NBT(:,k)/norm( NBT(:,k),2)))
    newlabels{i} = " ";
else
    newlabels{i} = labels(i);
end
end

labelnode(NBTplot,1:size(A,2), string(newlabels))



%Katz Graph

figure 


 Katzplot = plot(graph(A),  'Layout',"force", 'NodeCData', Katz(:,k))
 Katzplot.MarkerSize =  50*( Katz(:,k)/norm( Katz(:,k),2))
 Katzplot.NodeFontSize = 12
colorbar

for i=1:size(A,2)
if  Katz(i,k)  < 0.5*max( ( Katz(:,k)/norm( Katz(:,k),2)))
    newlabels{i} = " ";
else
    newlabels{i} = labels(i);
end
end

labelnode(Katzplot,1:size(A,2), string(newlabels))

