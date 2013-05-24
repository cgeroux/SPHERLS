#include "exception2.h"

const char* exception2::what()const throw(){
  return sMsg.c_str();
}
std::string exception2::getMsg(){
  return sMsg;
}
void exception2::setMsg(std::string sMsgIn){
  sMsg=sMsgIn;
}
void exception2::setCode(int nCodeIn){
  nCode=nCodeIn;
}
int exception2::getCode(){
  return nCode;
}
exception2::~exception2() throw(){
}
exception2::exception2():exception(){
}
exception2::exception2(std::string sMsg,int nCode):exception(),sMsg(sMsg),nCode(nCode){
}
exception2::exception2(std::string sMsg):exception(),sMsg(sMsg){
}
exception2::exception2(const exception2 &exception2In):exception(){//copy constructor
  setMsg(exception2In.sMsg);
  setCode(exception2In.nCode);
}
