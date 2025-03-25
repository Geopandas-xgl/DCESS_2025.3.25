function BorderLine(obj,varargin)
limX = get(obj, 'Xlim');
limY = get(obj, 'Ylim');
minX = limX(1);
maxX = limX(2);
minY = limY(1);
maxY = limY(2);
if length(varargin)>=1
line([maxX,maxX], [minY,maxY], "Color", "k","LineWidth",varargin{1});
line([minX,maxX], [maxY,maxY], "Color", "k","LineWidth",varargin{1});
else
line([maxX,maxX], [minY,maxY], "Color", "k","LineWidth",1);
line([minX,maxX], [maxY,maxY], "Color", "k","LineWidth",1);    
end
end