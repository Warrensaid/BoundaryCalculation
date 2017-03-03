tImage=zeros(xmax,ymax);
for n=1:xmax
    for m=1:ymax
        if mImage(n,m)==8
            tImage(n,m)=50;
        elseif mImage(n,m)==10
            tImage(n,m)=100;
        end
    end
end
figure
imagesc(tImage)