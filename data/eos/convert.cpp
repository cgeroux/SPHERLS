#include <fstream>
#include <iostream>
#include <cmath>
#include "eos.h"

int main(){
  //read in bob's eos table
  eos eosConvert;
  
  try{
    //eosConvert.readBobsAscii("ct_opac");
    //eosConvert.readAscii("Chris_new_EOS_opac");
    //eosConvert.writeAscii("eos2.txt");
    //eosConvert.writeBin("eos2");
    //eosConvert.writeBin("eosCTOPAC");
    //eosConvert.readAscii("eosCTOPAC_fixed.txt");
    //eosConvert.readAscii("eosY240Z002.txt");
    //eosConvert.readBin("eosY240Z002");
    eosConvert.readBin("eosCTOPAC");
    eosConvert.writeAscii("eosCTOPAC.txt");
    eosConvert.readBin("eosCTOPAC_old");
    eosConvert.writeAscii("eosCTOPAC_old.txt");
    eosConvert.readBin("eosNewY240Z002");
    eosConvert.writeAscii("eosNewY240Z002.txt");
    eosConvert.readBin("eosY240Z002");
    eosConvert.writeAscii("eosY240Z002.txt");
    eosConvert.readBin("eosY300Z002");
    eosConvert.writeAscii("eosY300Z002.txt");
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
