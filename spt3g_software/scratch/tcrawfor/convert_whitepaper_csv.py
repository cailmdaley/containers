import csv, glob

file1=glob.glob('*.csv')
authors1 = []
authors2 = []
titles = []

with open(file1[0]) as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    for row in reader:
        if row[0] == '1' or row[0] == '2':
            authors1.append(row[3])
            authors2.append(row[4])
            titles.append(row[6])

f2 = open(file1[0].replace('csv','bib'),'w')

for i in np.arange(len(authors1)):
    f2.write('@article{swp-'+(authors2[i].lower()).replace(' ','')+',\n')
    f2.write('author = "'+authors1[i]+' '+authors2[i]+' and others",\n')
    f2.write('title = "'+titles[i]+'",\n')
    f2.write('year = 2019,\n')
    f2.write('eprint = "Decadal Survey Science White Paper"\n')
    f2.write('}\n')
    f2.write('\n')

f2.close()
