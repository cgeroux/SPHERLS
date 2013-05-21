/** @file
  
  Header file for \ref exceptoin2.cpp
*/

#ifndef EXCEPTION2_H
#define EXCEPTION2_H

#include <exception>
#include <string>

//common error types
#define SYNTAX 0
#define INPUT 1
#define CALCULATION 2
#define OUTPUT 3
#define UNKNOWN 4

class exception2: public std::exception{
  private:
    //char* cMsg;
    std::string sMsg;
    int nCode;
  public:
    exception2();
    exception2(std::string sMsg,int nCode);
    exception2(std::string sMsg);
    exception2(const exception2 &exception2In);//copy constructor
    std::string getMsg();
    void setMsg(std::string sMsg);
    void setCode(int nCodeIn);
    int getCode();
    virtual ~exception2() throw();
};/**@class exception2
  Adds custom exception handling class
*/

#endif
