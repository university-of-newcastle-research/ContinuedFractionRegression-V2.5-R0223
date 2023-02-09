from sklearn.model_selection import train_test_split
import numpy as np
import os
import sys

if __name__ == '__main__':
    inputFilepath = sys.argv[1] #read input file path from the command-line

    filename_w_ext = os.path.basename(inputFilepath)
    fileName, file_extension = os.path.splitext(filename_w_ext)

    data=[]
    with open(inputFilepath, "r") as my_file:
        next(my_file)
        for line in my_file:
            list = line.split()
            data.append(list)
    my_file.close()

    npdata = np.array(data)  # convert array to numpy type array
    x_train, x_test = train_test_split(npdata, train_size=0.75, random_state=1)  # test_size=0.5(whole_data)

    nsamp, nfeat = x_train.shape
    headerTxt = str(nsamp)+"\t"+str(nfeat-1)
    outName = fileName+"-trn.in"
    np.savetxt(outName, x_train, header=headerTxt, comments='', delimiter='\t', fmt="%s")  # use exponential notation

    nSamp, nFeat = x_test.shape
    headerTxt = str(nSamp) + "\t" + str(nFeat-1)
    outName = fileName + "-Tst.in"
    np.savetxt(outName, x_test, header=headerTxt, comments='', delimiter='\t', fmt="%s")  # use exponential notation

    print(str(nsamp) + ", " + str(nfeat))
    print('Done!')
