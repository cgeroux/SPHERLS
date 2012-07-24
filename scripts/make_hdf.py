#!/usr/bin/env python
import sys
import paths
sys.path.append(paths.srcPath+"/pythonextensions/lib/python/")
import hdf
import dump

def main():
  dumpFile=dump.dump("./T5700_5x5_t00000000")
  dumpFile.printHeader()
  dumpFile.printVar(5)
if __name__ == "__main__":
  main()
'''
dims=[2,2,2]
data=[
[[1.0,1.1], [2.0,2.1]],
[[1.0,1.1], [2.0,2.1]]
]
hdf.open("test.hdf")
hdf.openData("test",hdf.DFNT_FLOAT64,dims)
hdf.writeData(data)
hdf.closeData()
hdf.openData("test2",hdf.DFNT_FLOAT64,dims)
hdf.closeData()
hdf.close()
'''