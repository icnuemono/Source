function [eGenoverA] = tempCalc(x,y) 
%constants
H=1;
SigmaX =4.5;              %x distribution
SigmaY =7.5;              %y distribution
xHat = [2.5:10:72.5];     %x locations
yHat = [2.75:16.5:52.25]; %y locations
sumAll = 0;
  for c=1:8
    for r=1:4
      %for each core sum the total heat contribution 
        sumAll = sumAll + exp(-((x-xHat(c)).^2)./(2*SigmaX^2)).*...
          exp(-((y-yHat(r)).^2)./(2*SigmaY^2)); 
    end
  end
  eGenoverA = sumAll;
end
