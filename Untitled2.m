a=[];
for i = 1:5
try
            GMModels = fitgmdist([1,2,3,4,5,6,7]',i,'Replicates',5);
            a(:,i)=cluster(GMModels,[1,2,3,4,5,6,7]');
        catch
            return
        end
end