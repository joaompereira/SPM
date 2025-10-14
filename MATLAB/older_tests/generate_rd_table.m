load results/testrd_table.mat

tablat = "Tensor,m,d,r,Extract (s),Power Method (s),Deflate (s),";
tablat = tablat + "Total time (s),Avg iter,Restarts,Error";

for i=1:nrows
    
    tablat = tablat + sprintf('T%d,%d,%d,%d,',...
            i, stat(i).n, stat(i).L, stat(i).R);
    for j=4:7
        tablat = tablat + sprintf('%5.2f,', stat(i).(fields{j}));
    end
    tablat = tablat + sprintf('%.0f,%d,', stat(i).(fields{8}),...
                                   stat(i).(fields{9}));
    tablat = tablat + sprintf('%1.2E\n', stat(i).(fields{10}));   
end

fileID = fopen('results/table.csv','w');
fprintf(fileID, "%s", tablat);
fclose(fileID);