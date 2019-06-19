
mins = 33;
secs = 40;


%jouletotal = 17819; laps = 3; %practice 1
jouletotal = 40596; laps = 7; %practice 2

distance = laps * 6388.10 / 3.28084;
secs = secs + mins * 60;

avgSpeed = distance ./ secs;

kwh = jouletotal / 3.6e6;
kmkwh = distance / 1000 / kwh;
mikwh = kmkwh / 1.609;