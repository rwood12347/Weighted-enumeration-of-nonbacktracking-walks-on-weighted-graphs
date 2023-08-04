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
load('D:\Fauci paper redo\Temporal Network\Matrices\Temporal_M.mat');
load('D:\Fauci paper redo\Temporal Network\Matrices\Temporal_bigR.mat');
load('D:\Fauci paper redo\Temporal Network\Matrices\Temporal_bigL.mat');
load('D:\Fauci paper redo\Temporal Network\Matrices\Temporal_As.mat');
load('D:\Fauci paper redo\Temporal Network\Matrices\Temporal_sqrtS.mat');
load('D:\Fauci paper redo\Temporal Network\Matrices\Temporal_rhos_KatzNBT.mat');
n = size(As{1},1);
m = size(M,1);
percent_max = 0.95

Katz_t = percent_max/abs(both_rhos(1));
NBT_t = percent_max/abs(both_rhos(2));

pencil = eye(size(M)) - NBT_t*M;
%resolvent = pencil \ eye(size(M));
resolvent = inv(pencil);
x_NBT = (eye(n) +NBT_t*bigL.'*sqrtS*resolvent*sqrtS*bigR)*ones(77,1);
    
    
%The following loop calculates Katz centrality with attentuation factor t
%for the time evolving network.
x_katz = ones(n,1);
for j = 0:99
    i = 100-j;
    x_katz = Katz_method(As{i}, Katz_t, x_katz);
end
%normalise the centrality vector

x_katz = x_katz / norm(x_katz,2);
x_NBT = x_NBT/norm(x_NBT,2);


[~,NBT_top_indices] = sort(x_NBT, 'descend');
NBT_top_indices = NBT_top_indices(1:10);

[~,Katz_top_indices] = sort(x_katz, 'descend');
Katz_top_indices = Katz_top_indices(1:10);

subplot(1,2,2) 

union_indices = union(Katz_top_indices,NBT_top_indices)
union_values = [ x_NBT(union_indices),x_katz(union_indices)];
union_labels = labels(union_indices);

barchart = bar(categorical(union_labels), union_values, 'LineStyle', "--", 'LineWidth', 2, 'BarWidth', 1)
barchart(1).LineStyle = ":";
legend('NBT Katz','Katz')
ax = gca;
% ax.FontSize = 36;
% ax.Legend.FontSize = 28;
ax.Legend.LineWidth = 1.5;
ax.Title.String = sprintf("Katz centrality values with attenuation factor t = %.2g/ρ", percent_max);
% ax.Title.FontSize = 36;
ax.YAxis.Label.String = 'Normalized centrality value';
% ax.XAxis.Label.String = 'Union of 10 most highly ranked nodes for Katz and NBT Katz';
hold on
