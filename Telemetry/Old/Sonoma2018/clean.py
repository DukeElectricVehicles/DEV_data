with open('cleaned.txt', 'w') as outfile:
    with open("sonoma2018day1raw.csv", "r") as ins:
        for i, line in enumerate(ins):
            line = line.strip()

            parts = line.split()
            parts = parts[:10]

            for p in parts:
                blah = float(p)
                if len(p) < 1:
                    print 'not enough chars', i

            if float(parts[3]) > 15:
                print "bad velo", i, parts[3], line
                continue

            new = ' '.join(parts)

            outfile.write(new + '\n')
