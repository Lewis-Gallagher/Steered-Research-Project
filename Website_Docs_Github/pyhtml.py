infile = "rn4_altpipeline.txt"
outfile = "cleaning5.txt"

delete_list = ["class=", "tg-yw4l"]
fin = open(infile)
fout = open(outfile, "w+")
for line in fin:
    for word in delete_list:
        line = line.replace(word,' ' )
    fout.write(line)
fin.close()
fout.close()
