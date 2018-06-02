function [same_pos,New_pow,adjacentAssemblies]=readQ_sym(powerDetectorFiles,assemblyPowerThreshold)
Power = readQ(powerDetectorFiles);
lengthQ_original = length(Power);
[Power, map] = formMap(Power, assemblyPowerThreshold);
nass=length(Power);
adjacentAssemblies = findAdjacentAssemblies(nass, map, lengthQ_original);
Sixths_index;
for m=1:size(Power,2)
Q=Power(:,m);
for i=1:6
    for j=1:length(sixths{i})
        for k=1:length(sixths{i}{j})
            same_pos{j}(i,k)=sixths{i}{j}(k);
            Q_same_pos{j}(i,k)=Q(indices(sixths{i}{j}(k),1));
        end
    end
end

for i=1:length(same_pos)
    for j=1:size(same_pos{i},2)
        Q_new_sixth{i}(j,1)=mean(Q_same_pos{i}(:,j));
        err{i}(j,1)=std(Q_same_pos{i}(:,j));
    end
end
        
%%
new_Q_rings=Q_new_sixth{1}(1);
for i=2:length(Q_new_sixth)
%     for j=1:length(Q_new_sixth{i})
        new_Q_rings=[new_Q_rings;repmat(Q_new_sixth{i},6,1)];
%     end
end
%%

for i=1:length(new_Q_rings)
    new_Q(indices(i,1),1)=new_Q_rings(i);
end

New_pow(:,m)=new_Q;
end
same_pos(1)=[];
same_pos=cell2mat(same_pos);
for i=1:size(same_pos,1)
    for j=1:size(same_pos,2)
        same_pos(i,j)=indices(same_pos(i,j),1);
    end
end

% plot_vect(adjacentAssemblies,new_Q/1e6)