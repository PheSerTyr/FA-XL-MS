database=open("BIOGRID-ORGANISM-Mus_musculus-4.4.233.tab3.txt","r")
database.readline()

xpair=[]

for line in database.readlines():
    words=line.split("\t")
    p=[words[23],words[26]]
    p.sort()
    xpair.append(tuple(p))

database.close()

xlfile=open("D:\\MS_DATA\\plus_20251218\\Testis_GluC_IAPI\\rep2\\1\\XL_Interlinks.tsv_q.tsv","r")
sw=open("D:\\MS_DATA\\plus_20251218\\Testis_GluC_IAPI\\rep2\\1\\XL_Interlinks.tsv_checked_IL.tsv","w")
sw.write(xlfile.readline().strip("\n")+"\tInteraction Check\n")

for line in xlfile.readlines():
    words=line.split("\t")
    p1=words[9].split("|")
    p2=words[22].split("|")
    interpair=[]

    intraflag=False

    for i in p1:
        if intraflag:
            break
        for j in p2:
            if i==j:
                sw.write(line.strip("\n")+"\tIntra(?)\n")
                intraflag=True
                break
            p=[i,j]
            p.sort()
            if tuple(p) in xpair:
                interpair.append("("+" & ".join(p)+")")
    
    if intraflag:
        continue
    if len(interpair)!=0:
        sw.write(line.strip("\n")+"\tYes\t"+";".join(interpair)+"\n")
    else:
        sw.write(line.strip("\n")+"\tNo\n")

sw.close()
xlfile.close()


print("end")

