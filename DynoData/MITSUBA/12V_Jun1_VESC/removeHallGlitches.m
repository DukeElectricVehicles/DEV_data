function [rpm] = removeHallGlitches(rpm)

    for i = 1:length(rpm) - 2%fix glitches in rpm readout
       if (rpm(i) > 0) && (rpm(i+2) > 0) && (rpm(i+1) == 0)
           rpm(i+1) = rpm(i);
       end
    end
    
    time = 1./rpm;
    
    for i = length(time):-1:10
        try
            if (time(i-1) < time(i)*.6)
                time(i-2:i-1) = time(i-2:i-1)*2;
            end
        catch e
            time(i-3:i+3)
        end
    end
    
    rpm = 1./time;
end