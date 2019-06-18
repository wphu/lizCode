import numpy as np
import h5py as h5
import glob
import os

def collect(groupName, dataSetName = None, itime=None, path="data",prefix="data"):
    file_list = glob.glob(os.path.join(path, prefix+"*.h5"))
    file_list.sort()

    if itime == None:
        i = 0
        for fileName in file_list:
            f = h5.File(fileName)
            if dataSetName == None:
                val = val = f[ groupName ]
            else:
                val = f[ groupName + dataSetName ]
            val = val[:]
            val = val[np.newaxis, :]
            if i == 0:
                data = val
            else:
                data = np.append(data, val, axis = 0)
            i = i + 1
            f.close()
        return data
    else:
            f = h5.File(file_list[itime])
            if dataSetName == None:
                val = val = f[ groupName ]
            else:
                val = f[ groupName + dataSetName ]
            val = val[:]
            f.close()
            return val   




#pn = collect("/Diagnostic", "particleNumber")     
