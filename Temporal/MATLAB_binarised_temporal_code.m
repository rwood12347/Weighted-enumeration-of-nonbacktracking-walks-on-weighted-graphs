clear all
labels = ["Fauci, Anthony";
"Haskins, Melinda";
"Selgrade, Sara";
"Crawford, Chase";
"Conrad, Patricia";
"Birx, Deborah";
"Tabak, Lawrence";
"Park, Alice";
"Emanuel, Ezekiel J.";
"Cassetti, Cristina";
"Routh, Jennifer";
"Folkers, Greg";
"Farrar, Jeremey";
"Frieden, Thomas";
"Awwad, David";
"Simonson, Stewart";
"Gao, George";
"Graham, Barney";
"Billet, Courtney";
"Lane, Cliff";
"Marston, Hilary";
"Dzau, Victor";
"Lerner, Andrea";
"Corey, Larry";
"Lapook, Jon";
"Jeffrey V Ravetch";
"Burklow, John";
"Collins, Francis";
"Anderson, Jennifer";
"Haynes, Barton";
"Stover, Kathy";
"Holland, Steven";
"Myles, Renate";
"Redfield, Robert";
"Giroir, Brett";
"Kadlec, Robert";
"Johnson, Martin";
"Auchincloss, Hugh";
"Harper, Jill";
"quinn, Thomas";
"Schneider, Johanna";
"Morens, David";
"Lipkin, Ian";
"Tromberg, Bruce";
"Denny, Joshua";
"Freire, Maria";
"Hatchett, Richard";
"Berkowitz, Avrahm";
"Tobias, Janet";
"Marks, Peter";
"Hahn, Stephen";
"Halkovich, Connie";
"Gallo, Robert";
"Fine, Amanda";
"Sabeti, Pardis";
"Miller, Katie";
"Shapiro, Neil";
"O'malley, Devin";
"Bright, Rick";
"Del rio, Carlos";
"Bertuzzi, Stefano";
"Glass, Roger";
"Kilmarx, Peter";
"Berkley, Seth";
"Mathewson, Herbert";
"Holdren, John P";
"Dirks, John";
"Lewis M Drusin";
"Eigabalawy, Nadia";
"Kelsall, Brian";
"Slavkin, Harold";
"Nader, Ralph";
"Galvani, Alison";
"Feinberg, Mark";
"Ghafoor, Azam";
"Emini, Emilio";
"Mundel, Trevor"
];
load('D:\Fauci paper redo\Temporal Network\Matrices\ binarised Temporal_M.mat');
load('D:\Fauci paper redo\Temporal Network\Matrices\ binarised Temporal_bigR.mat');
load('D:\Fauci paper redo\Temporal Network\Matrices\ binarised Temporal_bigL.mat');
load('D:\Fauci paper redo\Temporal Network\Matrices\ binarised Temporal_As.mat');
load('D:\Fauci paper redo\Temporal Network\Matrices\ binarised Temporal_sqrtS.mat');
load('D:\Fauci paper redo\Temporal Network\Matrices\ binarised Temporal_rhos_KatzNBT.mat');

n = size(As{1},1);
m = size(M,1);

%Make the ranges of attenuation factors
no_steps = 100;
Katz_upper_limit = 0.99/abs(both_rhos(1));
NBT_upper_limit = 0.99/abs(eigs(M,1));

Katz_t_range = 0:Katz_upper_limit/no_steps:Katz_upper_limit;
NBT_t_range = 0:NBT_upper_limit/no_steps:NBT_upper_limit;




%loop over the above ranges
for k = 1:no_steps+1
    Katz_t = Katz_t_range(k);
    NBT_t = NBT_t_range(k);
    resolvent = inv(speye(size(M)) - NBT_t*M);
    x_NBT = (speye(n) + NBT_t*bigL.'*resolvent*bigR)*ones(77,1);
    
    
%The following loop calculates Katz centrality with attentuation factor t
%for the time evolving network.
x_katz = ones(n,1);
resolvent = eye(n);
for j = 0:99
    resolvent = resolvent*inv(eye(n) - Katz_t*As{j+1});
end
%normalise the centrality vector
x_katz = resolvent*ones(n,1);

x_katz = x_katz / norm(x_katz,2);
x_NBT = x_NBT/norm(x_NBT,2);

    
Katz(:,k) = x_katz;
NBT(:,k) = x_NBT;
end


no_top_display = 10

% top indices for plot based on final centrality vector
[~,Katz_top_indices] = sort(Katz(:,no_steps+1), 'descend');
Katz_top_indices = Katz_top_indices(1:no_top_display);
[~,NBT_top_indices] = sort(NBT(:,no_steps+1), 'descend');

NBT_top_indices = NBT_top_indices(1:no_top_display);


%plot Katz centrality
marker_styles = ["o", "+", "*", ".","x","_","|", "square", "diamond", "^"];
line_styles = ["-", "--", ":", "-.","-", "--", ":", "-.","-", "--" ];

colours = {"#0072BD", "#D95319", "#EDB120", "#7E2F8E", "#77AC30", "#4DBEEE", "#A2142F","#0072BD", "#D95319", "#EDB120"};
for node = 1:no_top_display
  

 Katz_top_indices(node);
    XVec = Katz_t_range*abs(both_rhos(1));
    YVec = Katz(Katz_top_indices(node), :);
 if node == 6
     plot(XVec, YVec, 'LineWidth', 2, 'LineStyle', line_styles(node),'color', string(colours(node)), "Marker", marker_styles(node), 'MarkerIndices',10:10:length(XVec),  'MarkerFaceColor', string(colours(node)), 'MarkerSize',12)
 else
    
    if node == 4
        
    plot(XVec, YVec, 'LineWidth', 2, 'LineStyle', line_styles(node),'color', string(colours(node)) , "Marker", marker_styles(node), 'MarkerIndices',10:10:length(XVec),  'MarkerSize',36)
    else
       if node == 3
           plot(XVec, YVec, 'LineWidth', 2, 'LineStyle', line_styles(node),'color', string(colours(node)), "Marker", marker_styles(node), 'MarkerIndices',10:10:length(XVec),  'MarkerFaceColor', string(colours(node)),'MarkerSize',12)
       else
           
    plot(XVec, YVec, 'LineWidth', 2, 'LineStyle', line_styles(node), 'color', string(colours(node)) , "Marker", marker_styles(node), 'MarkerIndices',10:10:length(XVec),'MarkerSize',12)
       end
    end
 end
title('Katz Centrality for the Binarised Temporal Network')
ylabel('Normalized Katz centrality value');
xlabel('Attenuation factor as a percentage of 1/max_i (ρ(A^{[i]}))')

hold on
end
legend(labels(Katz_top_indices))



figure


I = eye(10);
P = [I(:,5), I(:,3), I(:,4), I(:,1), I(:,6), I(:,9), I(:,8), I(:,7), I(:,10), I(:,2)];
NBT_top_indices = P.'*NBT_top_indices;

for node = 1:no_top_display
    
    
    %
    
    XVec =NBT_t_range*abs(both_rhos(2));
    YVec =  NBT(NBT_top_indices(node), :);
 if node == 6
     plot(XVec, YVec, 'LineWidth', 2, 'LineStyle', line_styles(node),'color', string(colours(node)), "Marker", marker_styles(node), 'MarkerIndices',10:10:length(XVec),  'MarkerFaceColor', string(colours(node)), 'MarkerSize',12)
 else
    
    if node == 4
        
    plot(XVec, YVec, 'LineWidth', 2, 'LineStyle', line_styles(node),'color', string(colours(node)) , "Marker", marker_styles(node), 'MarkerIndices',10:10:length(XVec),  'MarkerSize',36)
    else
       if node == 3
           plot(XVec, YVec, 'LineWidth', 2, 'LineStyle', line_styles(node),'color', string(colours(node)), "Marker", marker_styles(node), 'MarkerIndices',10:10:length(XVec),  'MarkerFaceColor', string(colours(node)),'MarkerSize',12)
       else
           
    plot(XVec, YVec, 'LineWidth', 2, 'LineStyle', line_styles(node), 'color', string(colours(node)) , "Marker", marker_styles(node), 'MarkerIndices',10:10:length(XVec),'MarkerSize',12)
       end
    end
 end
    
    
    
    
    %
title('Non-backtracking Katz Centrality for the Binarised Temporal Network');
ylabel('Normalized NBT Katz centrality value');
xlabel('Attenuation factor as a percentage of 1/ρ(M)')
hold on
end


 legend(labels(NBT_top_indices));
 
 
 [labels(Katz_top_indices),labels(NBT_top_indices)];
 
 
 
 

 
 
 
