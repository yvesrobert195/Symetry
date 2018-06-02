[ind,r]=find_rings_and_numbers(adjacentAssemblies);

[~,ord]=sort(ind(:,2));
indices=ind(ord,:);
sixths=cell(6,1);
for i=1:6
    sixths{i}{1}=1;
end
ring=NaN;
i=2;
while i<length(indices)
    index=indices(i,1);
    if r(index)~=ring
        ring=r(index);
        for j=1:6
            sixths{j}{ring}=[];
        end
    end
    for j=1:6
        for k=1:ring-1
            if i>length(indices)
                break
            else
                sixths{j}{ring}(end+1)=ind(i,1);
                i=i+1;
            end
        end
    end
end

summ=1;
for i=1:ring-1
for j=1:(6*(i-1))
summ=summ+1;
end
end
last_ring=(nass-summ)/6;
n=1;
for i=1:6
    sixths{i}{ring}=[];
    for j=1:last_ring
        sixths{i}{ring}(j)=ind(summ+n,1);
        n=n+1;
    end
end


clear index i j k ring summ last_ring