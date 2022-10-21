import gzip
import sys

file0=gzip.open(sys.argv[1],'rb')

file1=gzip.open(sys.argv[2], 'wb')
file2=gzip.open(sys.argv[3], 'wb')

counter=0
for line in file0:
    if counter<4:
        file1.write(line)
        counter=counter+1
    else:
        if counter==7:
            file2.write(line)
            counter=0
        else:
            file2.write(line)
            counter=counter+1
file0.close()
file1.close()
file2.close()
