# helper script to convert football.graphwrong into football-ground-truth.graph

fin = open("football.graphwrong", 'r')
count = 0
fout = open("football-ground-truth.graph", 'w')

for line in fin:
	# get the line, remove newline
	l = line.replace("\n", "")

	# if not empty string
	if l:
		i = int(l)

		if i != count:
			sout = str(count) + " " + str(i) + "\n"
			print sout,
			fout.write(sout)
			count = count + 1

fin.close()
fout.close()