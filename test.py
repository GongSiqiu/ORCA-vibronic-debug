#!/usr/bin/python
# -*- coding: UTF-8 -*-
import numpy as np
from strawberryfields.apps import data, qchem
from strawberryfields.apps import plot

def ReadNM(hess_file):
    with open(hess_file, 'r') as inp:
        for line in inp:
            if not "$normal_modes" in line: # Flag for the Hessian group
                continue
            else:
                for data in inp:
                    ##--- Important Parameters
                    hess_size = int(data.split()[0])
                    n5_sets = hess_size//5 # Number of 5 columns sets.
                    rest = hess_size%5
                    if rest == 0:
                        n5_sets = n5_sets-1 # If the size is divisible by 5, then the sets would be larger.
                        rest = 5
                    ## Create the HESSIAN array
                    hess_data = np.zeros((hess_size,hess_size))
                    next(inp) # Jump the second line.
                    break
                #for data in inp: # Jump the second line (improve this)
                #    break
                count = 0
                while count <= n5_sets:
                    count2 = 0
                    ## First get the data from all sets of 5 columns except the last one if rest
                    # is not zero.
                    if count != n5_sets:
                        for data in inp:
                            columns = [ int(num) for num in range(count*5,(count)*5+5) ]
                            data = data.split()
                            n_line = int(data[0])
                            if  count2 >= hess_size:
                                break
                            for i in range(5):
                                j = columns[i]
                                hess_data[n_line][j] = float(data[i+1])
                            count2 += 1
                    ## This is the condition for the last sets of column that will not be zero
                    # if the size of the HESS matrix is not divisible by 5.
                    else:
                        for data in inp:
                            columns = [ int(num) for num in range(count*5,(count)*5+rest) ]
                            if data.strip() == '':
                                break
                            data = data.split()
                            n_line = int(data[0])
                            if  count2 >= hess_size:
                                break
                            for i in range(rest):
                                j = columns[i]
                                hess_data[n_line][j] = float(data[i+1])
                            count2 += 1
                    count += 1
    return hess_data

def ReadFre(hess_file):
    with open(hess_file, 'r') as inp:
        for line in inp:
            if not "$vibrational_frequencies" in line: # Flag for the Hessian group
                continue
            else:
                for data in inp:
                    ##--- Important Parameters
                    hess_size = int(data.strip())
                    break
            fre=np.zeros(hess_size)
            count=1
            for data in inp:
                #columns = [ int(num) for num in range(count*2,2) ]
                data = data.split()
                fre[count-1] = data[1]
                # for i in range(2):
                #     j = columns[i]
                #     hess_data[n_line][j] = float(data[i+1])
                count+=1
                if count>hess_size:
                    break
            return fre

def ReadAtoms(hess_file):
    with open(hess_file, 'r') as inp:

        for line in inp:
            if not "$atoms" in line: # Flag for the Hessian group
                continue
            else:
                for data in inp:
                    ##--- Important Parameters
                    hess_size = int(data.strip())
                    break
            atom=np.zeros((hess_size,4))
            count=0
            for data in inp:
                #columns = [ int(num) for num in range(count*2,2) ]
                data = data.split()
                for i in range(4):
                    atom[count][i]=data[i+1]
                count+=1
                if count>=hess_size:
                    break
            return atom


