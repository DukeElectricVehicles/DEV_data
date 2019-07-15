
mins = 33;
secs = 51;

jouletotal = 38264; laps = 7; %world record


distance = laps * 6388.10 / 3.28084;
secs = secs + mins * 60;

avgSpeed = distance ./ secs;

kwh = jouletotal / 3.6e6;
kmkwh = distance / 1000 / kwh;
mikwh = kmkwh / 1.60934;

wh100km = (kwh)/(distance/100e6);

%distance = 1;
%jouletotal = 2.9214;
mlgaseq = jouletotal / (0.7646 * 42.9e3); %ml gas eq
mpg = distance / (1609.34) / (mlgaseq / 3785.41);