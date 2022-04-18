function DrawDistanceTree(ps)

Y = pdist(ps);
z = linkage(Y);

h = dendrogram(z);
set(h,'LineWidth',2);
set(h,'Color','k');
xlabel('Neuron Index');
ylabel('Neuron Distance');
set(gca,'FontSize',18);

end