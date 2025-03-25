function ticklabels(axis,ticks,blank,retract,varargin)

% aixs: "x" or "y"
% ticks: ticks sequence
% blank: 1 or 2  number of interval between adjacent label;
% retract: "T" or "F" The label whether retract one blank? 
% varargin: Int e.g. "2" represent 2 decimal places  
    
    max1 = max(ticks);
    min1 = min(ticks);
    inter =  blank * (ticks(2) - ticks(1));
    if retract == "F"
        if length(varargin)>=1
        ary = round((min1:inter:max1),varargin{1});
        else
        ary = min1:inter:max1;
        end
        if blank == 1
            outtl = strings(size(ary, 1), size(ary, 2));
            for i = 1: size(ary,1)
            outtl(i,:) = ary(i,:);
            end
        elseif blank == 2
           outtl = strings(size(ary, 1), 2*size(ary, 2));
            for i = 1:size(ary, 1)
                outtl(i, 1:2:end-1) = ary(i, :);
                outtl(i, 2:2:end) = ' ';
            end
        end
     
    elseif retract == "T"
         min2 = min1+inter/blank;
          if length(varargin)>=1
            ary = round((min2:inter:max1),varargin{1});
         else
            ary = min2:inter:max1;
         end
         outtl = strings(size(ary, 1), 2*size(ary, 2));
            for i = 1:size(ary, 1)
                outtl(i, 1:2:end-1) = ' ';
                outtl(i, 2:2:end) = ary(i, :);
            end
    end
    
    tl = cellstr(outtl);
    if axis == "x"
    xticklabels(outtl);
    elseif axis == "y"
    yticklabels(outtl);
    end     
end