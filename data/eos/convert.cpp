#include <fstream>
#include <iostream>
#include <cmath>
#include "eos.h"

int main(){
  //read in bob's eos table
  try{
    eos eosConvert;
    //eosConvert.readBobsAscii("ct_opac");
    //eosConvert.readAscii("Chris_new_EOS_opac");
    //eosConvert.writeAscii("eos2.txt");
    //eosConvert.writeBin("eos2");
    //eosConvert.writeBin("eosCTOPAC");
    eosConvert.readAscii("eosCTOPAC_fixed.txt");
    eosConvert.writeBin("eosCTOPAC");
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
  //ofOut.open("pressure.dat");
  //ofOut.precision(16);
  //ofOut.unsetf(std::ios::fixed);
  //ofOut.setf(std::ios::scientific);
  //for(int i=0;i<eosConvert.nNumRho;i++){
  //  for(int j=0;j<eosConvert.nNumT;j++){
  //    ofOut<<eosConvert.dLogRhoMin+eosConvert.dLogRhoDelta*double(i)<<" "
  //      <<eosConvert.dLogTMin+eosConvert.dLogTDelta*double(j)<<" "
  //      <<eosConvert.dLogP[i][j]<<std::endl;
  //  }
  //  ofOut<<std::endl;
  //}
  //ofOut.close();
  
  ofOut.open("opacity.dat");
  ofOut.precision(16);
  ofOut.unsetf(std::ios::fixed);
  ofOut.setf(std::ios::scientific);
  for(int i=0;i<eosConvert.nNumRho;i++){
    for(int j=0;j<eosConvert.nNumT;j++){
      ofOut<<eosConvert.dLogRhoMin+eosConvert.dLogRhoDelta*double(i)<<" "
        <<eosConvert.dLogTMin+eosConvert.dLogTDelta*double(j)<<" "
        <<eosConvert.dLogKappa[i][j]<<std::endl;
    }
    ofOut<<std::endl;
  }
  ofOut.close();
  
  //ofOut.open("energy.dat");
  //ofOut.precision(16);
  //ofOut.unsetf(std::ios::fixed);
  //ofOut.setf(std::ios::scientific);
  //for(int i=0;i<eosConvert.nNumRho;i++){
  //  for(int j=0;j<eosConvert.nNumT;j++){
  //    ofOut<<eosConvert.dLogRhoMin+eosConvert.dLogRhoDelta*double(i)<<" "
  //      <<eosConvert.dLogTMin+eosConvert.dLogTDelta*double(j)<<" "
  //      <<eosConvert.dLogE[i][j]<<std::endl;
  //  }
  //  ofOut<<std::endl;
  //}
  //ofOut.close();*/
  return 0;
}
