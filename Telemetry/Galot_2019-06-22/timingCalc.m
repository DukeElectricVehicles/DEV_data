laps = 3;
distance = laps * 6388.10 / 3.28084;
v = 6.706;
startupTime = 15;
endTime = 20;

totTime = distance / v;
[totSec, totMin] = secToMin(totTime);

timeLap = (totTime - (startupTime + endTime)) / laps;

(totTime / laps) / timeLap  * 6.706

markIdx = 1:laps;
markTime = markIdx .* timeLap + startupTime;

for i = markIdx
   [lapSec, lapMin] = secToMin(markTime(i));
   fprintf('End lap %d time: %d min %.1f sec \n', i, lapMin, lapSec);
end

function [sec, min] = secToMin(secIn)
    min = floor(secIn/60);
    sec = secIn - min * 60;
end