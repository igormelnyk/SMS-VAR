
%p - process ID
%N - number of training data
%Nproc - total number of processes

%indFlights - indices of flights to process
%md - mode data
%sd - sensor data
%par - parameters
            
load(['sim/data_', num2str(p), '.mat']);

tic
for i=1:N/Nproc    
    jtree = JT(md{i}, sd{i}, par, 'sum');     
    
    jtree.propagateData(md{i}, sd{i}, par);
    
    res = jtree.getCopyWithAnomData();    
    
    save(['sim/res_',  num2str(indFlights(i)), '.mat'], '-v6', 'res');    
    
    fprintf('Data Sequence %d\n', i);
end
toc

disp('===== script evaluate sequences done =====');

exit

