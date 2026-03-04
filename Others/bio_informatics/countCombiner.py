import sys
mydict = {}
for file in sys.argv[1:]:
    for line in open(file, 'r'):
        key,value = line.strip().split('\t')
        if key in mydict:
            mydict[key] = mydict[key] + '\t' + value
        else:
            mydict[key] = value
for key,value in mydict.items():
    print(key + '\t' + value + '\n')