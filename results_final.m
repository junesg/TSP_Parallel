nodes_6 = [9.584 19.818, 34.242,  27.387 87.365 215.118];
nodes_12 = [12.411 16.628 36.882,  87.377 77.012 267.311];
nodes_24 = [23.169, 25.377, 47.864, 82.388 107.225 521.990];

size1 = [101 105 195  280 443 1002];


nodes_36 = [32.319 41.669, 57.463  89.923 114.202];
size2 = [101 105 195  280 443];


h=plot(size1, nodes_6, size1, nodes_12, size1, nodes_24, size2, nodes_36,size1, 1/5*size1,'--');
xlabel('Problem Size');
ylabel('Run time/sec');
title('Run time vs Problem size');
set(h, 'linewidth',3);
legend('6 nodes','12 ndoes', '24 nodes', '36 nodes','reference 1/5*problemSize');



Mat = [nodes_6; nodes_12; nodes_24; nodes_36 0];



quality = [0.979, 0.986,0.986,0.991;
    0.970 0.971 0.971 0.972
    0.931 0.955 0.956 0.956
    0.969 0.975 0.969 0.969
    0.958 0.958 0.968 0.968]
node_s = [6, 12,24,36];

h=plot(node_s, quality(1,:),node_s, quality(2,:),node_s, quality(3,:),node_s, quality(4,:),node_s, quality(5,:));
xlabel('Number of processors');
ylabel('Quality of solution as compared to TSPLIB');
title('Quality of solution vs number of processors');
set(h, 'linewidth',3);
legend('lin105','rat195','eil101','a280','pcb442');

