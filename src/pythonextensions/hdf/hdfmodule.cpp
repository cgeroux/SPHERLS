#include <Python.h>
#include <iostream>
#include "mfhdf.h"
#include <sstream>
#include <vector>

int32 nDataID=FAIL;
int32 nRank=0;
int32* nDims=NULL;
static PyObject* HDFError;
int32 nFileID=FAIL;

bool walkListAndWriteData(PyObject* list, int* nPosition,int nDepth){
  
  //check that list is really a list
  if(!PyList_Check(list)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": expecting a list for argument, check that you have nested the lists appropriately "
      <<"for the dimensions of the current varible";
    PyErr_SetString(HDFError,ssTemp.str().c_str());
    return false;
  }
  
  //check that the list is long enough
  if(PyList_Size(list)!=nDims[nDepth]){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": size of list given "<<PyList_Size(list)
      <<" doesn't match dimensions of varible previously set as "<<nDims[nDepth];
    PyErr_SetString(HDFError,ssTemp.str().c_str());
    return false;
  }
  
  if(nDepth==nRank-1){//we are at the bottom of the list, write out data
    
    //loop over list and write data
    double *dTemp= new double[nDims[nDepth]];
    nPosition[nDepth]=0;
    for(int i=0;i<nDims[nDepth];i++){
      
      PyObject* temp=PyList_GetItem(list,i);
      
      //check to make sure we really are at the bottom
      if(PyList_Check(temp)){
        std::stringstream ssTemp;
        ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
          <<": too many nested lists for data of rank "<<nRank;
        PyErr_SetString(HDFError,ssTemp.str().c_str());
        return false;
      }
      

      dTemp[i]=PyFloat_AsDouble(temp);
      
      //check to see if we might have an error some where
      if(dTemp[i]==-1){
        if(PyErr_Occurred()!=NULL){
          return false;
        }
      }
    }
    
    //write data to file
    int* nDimsTemp=new int[nRank];
    for(int j=0;j<nRank;j++){
      nDimsTemp[j]=1;
      if (j==nDepth){
        nDimsTemp[j]=nDims[nDepth];
      }
    }
    
    //DEBUG write out position
    /*std::cout<<"nPosition=(";
    for(int j=0;j<nRank-1;j++){
      std::cout<<nPosition[j]<<",";
    }
    std::cout<<nPosition[nRank-1]<<") ";
    
    //DEBUG write out position
    std::cout<<"nDimsTemp=(";
    for(int j=0;j<nRank-1;j++){
      std::cout<<nDimsTemp[j]<<",";
    }
    std::cout<<nDimsTemp[nRank-1]<<") ";
    
    //DEBUG write out position
    std::cout<<"dTemp=(";
    for(int j=0;j<nDims[nDepth];j++){
      std::cout<<dTemp[j]<<",";
    }
    std::cout<<dTemp[nDims[nDepth]-1]<<") \n\n";*/
    
    /**Warning this function call (SDWritedata) seems to leak memory, not sure why, I have tried a number of different methods to ensure that the buffer, dTemp is freed. I have also made sure that the data opened with SDcreate is closed with SDendaccess. Perhaps saving the data and doing one write would help, I am not sure. In any case doing only one write would very likely speed this up.
    
    I have moved this write out side the very bottom loop, and while it may not have removed the memory leak totaly, it seems to have drastically reduced it's severity.*/
    int32 nStat=SDwritedata(nDataID,nPosition,NULL,nDimsTemp,(VOIDP)dTemp);
    delete [] nDimsTemp;
    delete [] dTemp;
    if(nStat==FAIL){
      std::stringstream ssTemp;
      ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
        <<": failed writting data to file at data location (";
      for(int j=0;j<nRank-1;j++){
        ssTemp<<nPosition[j]<<",";
      }
      ssTemp<<nPosition[nRank-1]<<")";
      PyErr_SetString(HDFError,ssTemp.str().c_str());
      return false;
    }
  }
  else{//not at the bottom of the list, keep going
    //for each list in the current list call walkListAndWriteData
    
    for(int i=0;i<nDims[nDepth];i++){
      
      //get sub list
      PyObject* temp=PyList_GetItem(list,i);
      nPosition[nDepth]=i;
      if(!walkListAndWriteData(temp,nPosition,nDepth+1)){
        return false;
      }
    }
  }
  return true;
}

/** Opens an HDF file*/
static PyObject* open(PyObject* self, PyObject* args){
  const char* fileName;
  if(!PyArg_ParseTuple(args,"s",&fileName))
    return NULL;
  
  nFileID=SDstart(fileName,DFACC_CREATE);
  if(nFileID==FAIL){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": error opening the hdf file\""<<fileName<<"\"\n";
    PyErr_SetString(PyExc_IOError,ssTemp.str().c_str());
  }
  
  return Py_BuildValue("i",nFileID);
}
static PyObject* close(PyObject* self, PyObject* args){
  
  int32 nFileIDLoc=FAIL;

  if(!PyArg_ParseTuple(args,"|i",&nFileIDLoc)){
    return NULL;
  }
  
  //check to make sure we have a good file ID to use
  if(nFileIDLoc==FAIL){//if a file ID wasn't given
    if(nFileID==FAIL){//check to see if a file has been opened previously and use that file it has
      //nothing to be done
      Py_INCREF(Py_None);
      return Py_None;
    }
    else{
      nFileIDLoc=nFileID;
    }
  }
  
  int32 nStat=SDend(nFileID);
  if(nStat==FAIL){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": closing HDF file failed";
    PyErr_SetString(HDFError,ssTemp.str().c_str());
    return NULL;
  }
  nFileID=FAIL;
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject* openData(PyObject* self, PyObject* args){
  
  int32 nFileIDLoc=FAIL;
  int32 nDataIDLoc=FAIL;
  int32 nType;
  PyObject* listDims;
  
  //check to see if we have a data stream open
  if(nDims!=NULL){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": must close opened data before opening another one";
    PyErr_SetString(HDFError,ssTemp.str().c_str());
    return NULL;
  }
  
  char* cVarName=NULL;
  if(!PyArg_ParseTuple(args,"siO|i",&cVarName,&nType,&listDims,&nFileIDLoc)){
    return NULL;
  }
  
  //check to make sure we have a good file ID to use
  if(nFileIDLoc==FAIL){//if a file ID wasn't given
    if(nFileID==FAIL){//check to see if a file has been opened previously and use that file it has
      std::stringstream ssTemp;
      ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
        <<": no hdf file opened, open a file by calling \"hdf.open\" with the file name as an argument";
      PyErr_SetString(HDFError,ssTemp.str().c_str());
      return NULL;
    }
    else{
      nFileIDLoc=nFileID;
    }
  }
  
  //check that listDims is really a list
  if(!PyList_Check(listDims)){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": expecting a list for argument 3";
    PyErr_SetString(HDFError,ssTemp.str().c_str());
    return NULL;
  }
  
  //get rank and dimensions
  nRank=PyList_Size(listDims);
  nDims=new int32[nRank];
  for(int i=0;i<nRank;i++){
    
    PyObject* temp=PyList_GetItem(listDims,i);
    if(!PyInt_Check(temp)){
      std::stringstream ssTemp;
      ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
        <<": expecting a list of integers for argument 3";
      PyErr_SetString(HDFError,ssTemp.str().c_str());
      return NULL;
    }
    int32 nDimTemp=PyInt_AsLong(temp);
    nDims[i]=nDimTemp;
  }
  
  //create the HDF data
  nDataIDLoc=SDcreate(nFileIDLoc,cVarName,nType,nRank,nDims);
  if(nDataIDLoc==FAIL){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": creating a new HDF data set for variable \""<<cVarName<<"\" failed\n";
    PyErr_SetString(HDFError,ssTemp.str().c_str());
    return NULL;
  }
  nDataID=nDataIDLoc;
  return Py_BuildValue("i",nDataIDLoc);
}
static PyObject* writeData(PyObject* self, PyObject* args){
  PyObject* list;
  if(!PyArg_ParseTuple(args,"O",&list)){
    return NULL;
  }
  
  //create and initialize starting position
  int* nPosition=new int[nRank];
  for(int i=0;i<nRank;i++){
    nPosition[i]=0;
  }
  
  if(!walkListAndWriteData(list,nPosition,0)){
    delete [] nPosition;
    return NULL;
  }
  
  //nStat=SDwritedata(nDataID,nStart,NULL,nDims,(VOIDP)dTemp);
  delete [] nPosition;
  Py_INCREF(Py_None);
  return Py_None;
}
static PyObject* closeData(PyObject* self, PyObject* args){
  int32 nDataIDLoc=FAIL;
  
  //first lets delete nDims if not NULL
  if(nDims!=NULL){
    delete [] nDims;
  }
  nDims=NULL;
  nRank=0;
  
  //see if a dataID was given
  if(!PyArg_ParseTuple(args,"|i",&nDataIDLoc))
    return NULL;
  
  //check to make sure we have a good data ID to use
  if(nDataIDLoc==FAIL){//if a data ID wasn't given
    if(nDataID==FAIL){//check to see if a data has been opened previously and use that data it has
      //nothing to do
      Py_INCREF(Py_None);
      return Py_None;
    }
    else{
      nDataIDLoc=nDataID;
    }
  }
  int32 nStat=SDendaccess(nDataID);
  if(nStat==FAIL){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": error closing data";
    PyErr_SetString(PyExc_IOError,ssTemp.str().c_str());
    return NULL;
  }
  nDataID=FAIL;
  Py_INCREF(Py_None);
  return Py_None;
}
static PyMethodDef HDFMethods[] ={
     {"open", open, METH_VARARGS, "Opens an HDF file"},
     {"openData", openData, METH_VARARGS, "Opens access for writting data to hdf file"},
     {"writeData", writeData, METH_VARARGS, "Writes data to an hdf file"},
     {"closeData", closeData, METH_VARARGS, "Closes access for writting data to hdf file"},
     {"close", close, METH_VARARGS, "Closes hdf file"},
     {NULL, NULL, 0, NULL}
};
PyMODINIT_FUNC inithdf(void){
  PyObject *m;
  m=Py_InitModule("hdf", HDFMethods);
  if(m==NULL){
    return;
  }
  
  HDFError=PyErr_NewException("hdf.error",NULL,NULL);
  Py_INCREF(HDFError);
  PyModule_AddObject(m,"error",HDFError);
  
  //add some needed HDF constants
  PyModule_AddIntConstant(m, "DFNT_FLOAT64", DFNT_FLOAT64);
}
