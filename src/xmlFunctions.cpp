#include "xmlFunctions.h"
#include "exception2.h"
#include <string>
#include "xmlParser.h"
#include <sstream>

void getXMLValue(XMLNode xParent, std::string sChild,int nIndex,int &nValue){
  XMLNode xTemp=xParent.getChildNode(sChild.c_str(),nIndex);
  if(xTemp.isEmpty()){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": no \""<<sChild<<"\" node found under \""<<xParent.getName()<<"\" node\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  if(!xTemp.getText()){//if the node contains no text
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": \""<<sChild<<"\" node contains no value\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  nValue=atoi(xTemp.getText());
}
void getXMLValue(XMLNode xParent, std::string sChild,int nIndex,unsigned int &nValue){
  XMLNode xTemp=xParent.getChildNode(sChild.c_str(),nIndex);
  if(xTemp.isEmpty()){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": no \""<<sChild<<"\" node found under \""<<xParent.getName()<<"\" node\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  if(!xTemp.getText()){//if the node contains no text
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": \""<<sChild<<"\" node contains no value\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  nValue=atoi(xTemp.getText());
}
void getXMLValue(XMLNode xParent, std::string sChild,int nIndex, double &dValue){
  XMLNode xTemp=xParent.getChildNode(sChild.c_str(),nIndex);
  if(xTemp.isEmpty()){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": no \""<<sChild<<"\" node found under \""<<xParent.getName()<<"\" node\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  if(!xTemp.getText()){//if the node contains no text
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": \""<<sChild<<"\" node contains no value\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  dValue=atof(xTemp.getText());
}
void getXMLValue(XMLNode xParent, std::string sChild,int nIndex, std::string &sValue){
  XMLNode xTemp=xParent.getChildNode(sChild.c_str(),nIndex);
  if(xTemp.isEmpty()){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": no \""<<sChild<<"\" node found under \""<<xParent.getName()<<"\" node\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  if(!xTemp.getText()){//if the node contains no text
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": \""<<sChild<<"\" node contains no value\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  sValue=xTemp.getText();
}
void getXMLValue(XMLNode xParent, std::string sChild,int nIndex, bool &bValue){
  XMLNode xTemp=xParent.getChildNode(sChild.c_str(),nIndex);
  if(xTemp.isEmpty()){
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": no \""<<sChild<<"\" node found under \""<<xParent.getName()<<"\" node\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  if(!xTemp.getText()){//if the node contains no text
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": \""<<sChild<<"\" node contains no value\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  
  //convert to all upper case
  std::string sTemp=xTemp.getText();
  for(unsigned int i=0;i<sTemp.size();i++){
    sTemp[i]=toupper(sTemp[i]);
  }
  if(sTemp.compare("TRUE")==0){
    bValue=true;
  }
  else if (sTemp.compare("FALSE")==0){
    bValue=false;
  }
  else{//not a proper value
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": node \""<<sChild<<"\" has an incorrect boolean value \n";
    throw exception2(ssTemp.str(),INPUT);
  }
}
void getXMLAttribute(XMLNode xParent,std::string sAttribute,int &nValue){
  //check to see there is an an sAttribute
  if(!xParent.isAttributeSet(sAttribute.c_str())){//no attribute
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": no \""<<sAttribute<<"\" attribute set for node \""<<xParent.getName()<<"\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  if(!xParent.getAttribute(sAttribute.c_str())){//no value
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": no value for \""<<sAttribute<<"\" attribute set for node \""<<xParent.getName()<<"\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  nValue=atoi(xParent.getAttribute(sAttribute.c_str()));
}
void getXMLAttribute(XMLNode xParent, std::string sAttribute, std::string &sValue){
  //check to see there is an sAttribute
  if(!xParent.isAttributeSet(sAttribute.c_str())){//no attribute
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": no \""<<sAttribute<<"\" attribute set for node \""<<xParent.getName()<<"\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  if(!xParent.getAttribute(sAttribute.c_str())){//no value
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": no value for \""<<sAttribute<<"\" attribute set for node \""<<xParent.getName()<<"\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  sValue=xParent.getAttribute(sAttribute.c_str());
}
void getXMLAttribute(XMLNode xParent, std::string sAttribute, double &dValue){
  //check to see there is an sAttribute
  if(!xParent.isAttributeSet(sAttribute.c_str())){//no attribute
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": no \""<<sAttribute<<"\" attribute set for node \""<<xParent.getName()<<"\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  if(!xParent.getAttribute(sAttribute.c_str())){//no value
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<": no value for \""<<sAttribute<<"\" attribute set for node \""<<xParent.getName()<<"\"\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  dValue=atof(xParent.getAttribute(sAttribute.c_str()));
}
int getXMLValueNoThrow(XMLNode xParent, std::string sChild,int nIndex,int &nValue){
  XMLNode xTemp=xParent.getChildNode(sChild.c_str(),nIndex);
  if(xTemp.isEmpty()){
    return 0;
  }
  if(!xTemp.getText()){//if the node contains no text
    return 0;
  }
  nValue=atoi(xTemp.getText());
  
  return 1;
}
int getXMLValueNoThrow(XMLNode xParent, std::string sChild,int nIndex,unsigned int &nValue){
  XMLNode xTemp=xParent.getChildNode(sChild.c_str(),nIndex);
  if(xTemp.isEmpty()){
    return 0;
  }
  if(!xTemp.getText()){//if the node contains no text
    return 0;
  }
  nValue=atoi(xTemp.getText());
  
  return 1;
}
int getXMLValueNoThrow(XMLNode xParent, std::string sChild,int nIndex, double &dValue){
  XMLNode xTemp=xParent.getChildNode(sChild.c_str(),nIndex);
  if(xTemp.isEmpty()){
    return 0;
  }
  if(!xTemp.getText()){//if the node contains no text
    return 0;
  }
  dValue=atof(xTemp.getText());
  return 1;
}
int getXMLValueNoThrow(XMLNode xParent, std::string sChild,int nIndex, std::string &sValue){
  XMLNode xTemp=xParent.getChildNode(sChild.c_str(),nIndex);
  if(xTemp.isEmpty()){
    return 0;
  }
  if(!xTemp.getText()){//if the node contains no text
    return 0;
  }
  sValue=xTemp.getText();
  return 1;
}
int getXMLValueNoThrow(XMLNode xParent, std::string sChild,int nIndex, bool &bValue){
  XMLNode xTemp=xParent.getChildNode(sChild.c_str(),nIndex);
  if(xTemp.isEmpty()){
    return 0;
  }
  if(!xTemp.getText()){//if the node contains no text
    return 0;
  }
  
  //convert to all upper case
  std::string sTemp=xTemp.getText();
  for(unsigned int i=0;i<sTemp.size();i++){
    sTemp[i]=toupper(sTemp[i]);
  }
  if(sTemp.compare("TRUE")==0){
    bValue=true;
  }
  else if (sTemp.compare("FALSE")==0){
    bValue=false;
  }
  else{//not a proper value
    return 0;
  }
  return 1;
}
int getXMLAttributeNoThrow(XMLNode xParent,std::string sAttribute,int &nValue){
  //check to see there is an an sAttribute
  if(!xParent.isAttributeSet(sAttribute.c_str())){//no attribute
    return 0;
  }
  if(!xParent.getAttribute(sAttribute.c_str())){//no value
    return 0;
  }
  nValue=atoi(xParent.getAttribute(sAttribute.c_str()));
  return 1;
}
int getXMLAttributeNoThrow(XMLNode xParent, std::string sAttribute, std::string &sValue){
  //check to see there is an an sAttribute
  if(!xParent.isAttributeSet(sAttribute.c_str())){//no attribute
    return 0;
  }
  if(!xParent.getAttribute(sAttribute.c_str())){//no value
    return 0;
  }
  sValue=xParent.getAttribute(sAttribute.c_str());
  return 1;
}
int getXMLAttributeNoThrow(XMLNode xParent, std::string sAttribute, double &dValue){
  //check to see there is an an sAttribute
  if(!xParent.isAttributeSet(sAttribute.c_str())){//no attribute
    return 0;
  }
  if(!xParent.getAttribute(sAttribute.c_str())){//no value
    return 0;
  }
  dValue=atof(xParent.getAttribute(sAttribute.c_str()));
  return 1;
}
XMLNode getXMLNode(XMLNode xParent, std::string sChild,int nIndex){
  XMLNode xTemp=xParent.getChildNode(sChild.c_str(),nIndex);
  if(xTemp.isEmpty()){//if no node found
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": no \""<<sChild<<"\" node found under \""<<xParent.getName()<<"\" node\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  return xTemp;
}
XMLNode getXMLNodeNoThrow(XMLNode xParent, std::string sChild,int nIndex){
  XMLNode xTemp=xParent.getChildNode(sChild.c_str(),nIndex);
  return xTemp;
}
XMLNode openXMLFile(std::string sFileName, std::string sStartNode){
  XMLNode xTemp=XMLNode::openFileHelper(sFileName.c_str(),sStartNode.c_str());
  if(xTemp.isEmpty()){//check that we got a starting node
    std::stringstream ssTemp;
    ssTemp<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__
      <<": error parsing the file \""<<sFileName<<"\" possibly no global \""
      <<sStartNode<<"\" node\n";
    throw exception2(ssTemp.str(),INPUT);
  }
  return xTemp;
}
XMLNode openXMLFileNoThrow(std::string sFileName, std::string sStartNode){
  XMLNode xTemp=XMLNode::openFileHelper(sFileName.c_str(),sStartNode.c_str());
  return xTemp;
}
