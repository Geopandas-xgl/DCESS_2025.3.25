function text_norm(CoeX,CoeY,letter,FontSize) 
limX = get(gca, 'Xlim');
limY = get(gca, 'Ylim');
minX = limX(1);
maxX = limX(2);
minY = limY(1);
maxY = limY(2);
text(minX+(maxX-minX)*CoeX,minY+(maxY-minY)*CoeY,letter,"FontSize",FontSize,"FontName", "Times New Roman");
end