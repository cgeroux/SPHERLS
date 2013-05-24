# distutils: language = c++
# distutils: sources = eos_tmp.cpp exception2.cpp
from libcpp.string cimport string

cdef extern from "exception2.h":
  cdef cppclass exception2:
    string sMsg
    int nCode
    
    exception2(string)except +
    
    string getMsg()
    setMsg(string)
    const char* what()
cdef class Exception2:
  cdef exception2 *thisptr
  def __cinit__(self,msg):
    self.thisptr=new exception2(msg)
  def __dealloc__(self):
    del self.thisptr
  def getMsg(self):
    return self.thisptr.getMsg()
  def what(self):
    return self.thisptr.what()
cdef extern from "eos.h":
  cdef cppclass eos:
    eos() except +
    
    int nNumRho
    int nNumT
    double dXMassFrac
    double dYMassFrac
    double dLogRhoMin
    double dLogRhoDelta
    double dLogTMin
    double dLogTDelta
    double **dLogP
    double **dLogE
    double **dLogKappa
    
    void readAscii(string) except +
    void readBobsAscii(string) except +
    void writeAscii(string) except +
    void readBin(string) except +
    void writeBin(string) except +
    double dGetPressure(double, double) except +
    double dGetEnergy(double, double) except +
    double dGetOpacity(double, double) except +
    double dDRhoDP(double, double) except +
    double dSoundSpeed(double, double) except +
    void getEKappa(double, double, double&, double&) except +
    void getPEKappa(double, double, double&, double&, double&) except +
    void getPEKappaGamma(double, double, double&, double&, double&, double&) except +
    void getPEKappaGammaCp(double, double, double&, double&, double&, double&, double&) except +
    void getPKappaGamma(double, double, double&, double&, double&) except +
    void gamma1DelAdC_v(double, double, double&, double&, double&) except +
    void getPAndDRhoDP(double, double, double&, double&) except +
    void getEAndDTDE(double, double, double&, double&) except +
    void getDlnPDlnTDlnPDlnPDEDT(double, double, double&, double&, double&) except +
cdef class Eos:
  cdef eos *thisptr      # hold a C++ instance which we're wrapping
  def __cinit__(self):
    self.thisptr = new eos()
  def __dealloc__(self):
    del self.thisptr
  def readAscii(self, fileName):
    self.thisptr.readAscii(fileName)
  def readBobsAscii(self, fileName):
    self.thisptr.readBobsAscii(fileName)
  def writeAscii(self, fileName):
    self.thisptr.writeAscii(fileName)
  def readBin(self, fileName):
    self.thisptr.readBin(fileName)
  def writeBin(self, fileName):
    self.thisptr.writeBin(fileName)
  def getPressure(self, dTemp, dRho):
    return self.thisptr.dGetPressure(dTemp,dRho)
  def getEnergy(self, dTemp, dRho):
    return self.thisptr.dGetEnergy(dTemp,dRho)
  def getOpacity(self, dTemp, dRho):
    return self.thisptr.dGetOpacity(dTemp,dRho)
  def getDRhoDP(self, dTemp, dRho):
    return self.thisptr.dDRhoDP(dTemp,dRho)
  def getSoundSpeed(self, dTemp, dRho):
    return self.thisptr.dSoundSpeed(dTemp,dRho)
  def getEKappa(self, dTemp, dRho):
    e=0.0
    kappa=0.0
    self.thisptr.getEKappa(dTemp,dRho,e,kappa)
    return [e,kappa]
  def getPEKappa(self, dTemp, dRho):
    p=0.0
    e=0.0
    kappa=0.0
    self.thisptr.getPEKappa(dTemp,dRho,p,e,kappa)
    return [p,e,kappa]
  def getPEKappaGamma(self, dTemp, dRho):
    p=0.0
    e=0.0
    kappa=0.0
    gamma=0.0
    self.thisptr.getPEKappaGamma(dTemp,dRho,p,e,kappa,gamma)
    return [p,e,kappa,gamma]
  def getPEKappaGammaCp(self, dTemp, dRho):
    p=0.0
    e=0.0
    kappa=0.0
    gamma=0.0
    cp=0.0
    self.thisptr.getPEKappaGammaCp(dTemp,dRho,p,e,kappa,gamma,cp)
    return [p,e,kappa,gamma,cp]
  def getPKappaGamma(self, dTemp, dRho):
    p=0.0
    kappa=0.0
    gamma=0.0
    self.thisptr.getPKappaGamma(dTemp,dRho,p,kappa,gamma)
    return [p,kappa,gamma]
  def getGamma1DelAdC_v(self, dTemp, dRho):
    gamma1=0.0
    delAd=0.0
    cv=0.0
    self.thisptr.gamma1DelAdC_v(dTemp,dRho,gamma1,delAd,cv)
    return[gamma1,delAd,cv]
  def getPAndDRhoDP(self, dTemp, dRho):
    p=0.0
    DRhoDP=0.0
    self.thisptr.getPAndDRhoDP(dTemp,dRho,p,DRhoDP)
    return [p,DRhoDP]
  def getEAndDTDE(self, dTemp, dRho):
    e=0.0
    dTde=0.0
    self.thisptr.getEAndDTDE(dTemp,dRho,e,dTde)
    return [e,dTde]
  def getDlnPDlnTDlnPDlnPDEDT(self, dTemp, dRho):
    DlnPDlnT=0.0
    DlnPDlnRho=0.0
    DEDT=0.0
    self.thisptr.getDlnPDlnTDlnPDlnPDEDT(dTemp,dRho,DlnPDlnT,DlnPDlnRho,DEDT)
    return [DlnPDlnT,DlnPDlnRho,DEDT]
