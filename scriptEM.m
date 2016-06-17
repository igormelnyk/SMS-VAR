
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
    
    ll = jtree.run();
    logL = ll;
    
    jt = jtree.getPartialCopy(); 
    
    fprintf('Data Sequence %d\n', i);
    
    save(['sim/res_', num2str(indFlights(i)), '.mat'], '-v6', 'jt', 'logL');
end
toc

disp('====== script EM done ======');

exit