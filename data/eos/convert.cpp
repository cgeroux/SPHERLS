#include <fstream>
#include <iostream>
#include <cmath>
#include "eos.h"

int main(int argc, char* argv[]){
  //read in bob's eos table
  eos eosConvert;
  
  if (argc!=2){
    std::cout<<"expecting only one input file name"<<std::endl;
    return 0;
  }
  try{
    std::string sFileNameFrom=argv[1];
    std::string sExtension=sFileNameFrom.substr(sFileNameFrom.length()-4,4);
    if(sExtension!=".txt"){
      std::string sFileNameTo=sFileNameFrom+".txt";
      std::cout<<"converting binary eos file \""<<sFileNameFrom<<" to ascii file \""<<sFileNameTo
        <<"\"\n";
      eosConvert.readBin(sFileNameFrom);
      eosConvert.writeAscii(sFileNameTo);
    }
    else{
      std::string sFileNameTo=sFileNameFrom.substr(0,sFileNameFrom.length()-4);
      std::cout<<"converting ascii eos file \""<<sFileNameFrom<<" to binary file \""<<sFileNameTo
        <<"\"\n";
      eosConvert.readAscii(sFileNameFrom);
      eosConvert.writeBin(sFileNameTo);
    }
  }
  //error handeling
  catch(exception2& eTemp){
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<eTemp.getMsg();
  }
  catch(std::exception& eTemp){
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<"Standard exception:"<<eTemp.what()<<std::endl;
  }
  catch(...){
    std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<":"<<"main: unknown error\n";
  }
  /*
  std::ofstream ofOut;
  ofOut.open("pressure.dat");
  ofOut.precision(16);
  ofOut.unsetf(std::ios::fixed);
  ofOut.setf(std::ios::scientific);
  double dCurLogRho;
  double dCurLogT;
  //std::cout<<eosConvert.nNumRho<<std::endl;
  //std::cout<<eosConvert.nNumT<<std::endl;
  for(int i=0;i<eosConvert.nNumRho;i++){
    for(int j=0;j<eosConvert.nNumT;j++){
      
      dCurLogRho=eosConvert.dLogRhoMin+eosConvert.dLogRhoDelta*double(i);
      dCurLogT=eosConvert.dLogTMin+eosConvert.dLogTDelta*double(j);
      if(dCurLogRho>=-4.3&&dCurLogRho<=-3.89 &&dCurLogT>=5.66&&dCurLogT<=5.676){
        ofOut<<eosConvert.dLogP[i][j]<<" ";
      }
    }
    if(dCurLogRho>=-4.3&&dCurLogRho<=-3.89){
      ofOut<<std::endl;
    }
  }
  ofOut.close();
  
  ofOut.open("opacity.dat");
  ofOut.precision(16);
  ofOut.unsetf(std::ios::fixed);
  ofOut.setf(std::ios::scientific);
  for(int i=0;i<eosConvert.nNumRho;i++){
    for(int j=0;j<eosConvert.nNumT;j++){
      dCurLogRho=eosConvert.dLogRhoMin+eosConvert.dLogRhoDelta*double(i);
      dCurLogT=eosConvert.dLogTMin+eosConvert.dLogTDelta*double(j);
      if(dCurLogRho>=-4.3&&dCurLogRho<=-3.89 &&dCurLogT>=5.66&&dCurLogT<=5.676){
        ofOut<<eosConvert.dLogKappa[i][j]<<" ";
      }
    }
    if(dCurLogRho>=-4.3&&dCurLogRho<=-3.89){
      ofOut<<std::endl;
    }
  }
  ofOut.close();
  
  ofOut.open("energy.dat");
  ofOut.precision(16);
  ofOut.unsetf(std::ios::fixed);
  ofOut.setf(std::ios::scientific);
  for(int i=0;i<eosConvert.nNumRho;i++){
    for(int j=0;j<eosConvert.nNumT;j++){
      dCurLogRho=eosConvert.dLogRhoMin+eosConvert.dLogRhoDelta*double(i);
      dCurLogT=eosConvert.dLogTMin+eosConvert.dLogTDelta*double(j);
      //std::cout<<"("<<dCurLogRho<<","<<dCurLogT<<")"<<std::endl;
      if(dCurLogRho>=-4.3&&dCurLogRho<=-3.89 &&dCurLogT>=5.66&&dCurLogT<=5.676){
        std::cout<<"("<<dCurLogRho<<","<<dCurLogT<<") ";
        ofOut<<eosConvert.dLogE[i][j]<<" ";
      }
    }
    if(dCurLogRho>=-4.3&&dCurLogRho<=-3.89){
      ofOut<<std::endl;
      std::cout<<std::endl;
    }
  }
  ofOut.close();*/
  return 0;
}
