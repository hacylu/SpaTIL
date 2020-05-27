function [ feat,featDescription ] = getIntersectedArea( lympMembers,lympNodes,nonLympMembers,nonLympNodes )
%EXTRACTINTERSECTEDAREA Summary of this function goes here
%   Detailed explanation goes here
n=3;

[lympClusterCent,lymPolygons]=getClusterCentroids(lympMembers,lympNodes);
[nonLympClusterCent,nonLymPolygons]=getClusterCentroids(nonLympMembers,nonLympNodes);

cent=[lympClusterCent;nonLympClusterCent];
polyg=[lymPolygons nonLymPolygons];
d=squareform(pdist(cent));
d(d==0)=Inf;

numClust=length(cent);
intersectedArea=0;

if ~isempty(d)
    for i=1:numClust
        for j=1:n
            [~,ind]=min(d(i,:));
            d(i,ind)=Inf;
            
            polyarray1 = polyshape(polyg{i});
            polyarray2 = polyshape(polyg{ind});
            polyout = intersect([polyarray1 polyarray2]);
            intersectedArea=intersectedArea+polyarea(polyout.Vertices(:,1),polyout.Vertices(:,2));
            
            %[xp1,yp1]=poly2cw(pol1(:,1),pol1(:,2));
            %[xp2,yp2]=poly2cw(pol2(:,1),pol2(:,2));
            %[xi, yi] = polybool('intersection',xp1,yp1,xp2,yp2);
            %intersectedArea=intersectedArea+polyarea(xi,yi);
        end
    end
end

feat=intersectedArea/2;
featDescription={'intersectedArea'};

end