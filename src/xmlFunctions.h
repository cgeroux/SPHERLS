#ifndef XMLFUNCTIONS_H
#define XMLFUNCTIONS_H

#include <string>
#include "xmlParser.h"

void getXMLValue(XMLNode xParent, std::string sChild,int nIndex,int &nValue);
void getXMLValue(XMLNode xParent, std::string sChild,int nIndex, unsigned int &nValue);
void getXMLValue(XMLNode xParent, std::string sChild,int nIndex, double &dValue);
void getXMLValue(XMLNode xParent, std::string sChild,int nIndex, std::string &sValue);
void getXMLValue(XMLNode xParent, std::string sChild,int nIndex, bool &bValue);
void getXMLAttribute(XMLNode xParent, std::string sAttribute, int &nValue);
void getXMLAttribute(XMLNode xParent, std::string sAttribute, std::string &sValue);/**<
  Attempts to get the attribute named by \c sAttribute in the parent node \c xParent and returns its
  value in \c sValue. If no attribute is found it throws an exception letting the using know that 
  none was found in the node it was trying to check.
  */
void getXMLAttribute(XMLNode xParent, std::string sAttribute, double &dValue);/**<
  Attempts to get the attribute named by \c sAttribute in the parent node \c xParent and returns its
  value in \c dValue. If no attribute is found it throws an exception letting the using know that 
  none was found in the node it was trying to check.
  */
int getXMLValueNoThrow(XMLNode xParent, std::string sChild,int nIndex,int &nValue);/**<
  A version of \ref getXMLValue for integers that doesn't throw an exception if no xml node is 
  found. It returns 0 if node is either empty or not found, otherwise it returns 1 to indicate
  that the node exists and that it is not empty.
  
  @param[in] xParent the parent \ref XMLNode containing the node which contains the integer value
    being seeked.
  @param[in] sChild the name of the \ref XMLNode containing the integer value being seeked.
  @param[in] nIndex the index of the child node i.e. if there is more than one node of the same
    name and the second node is desired then nIndex should be 1 and nIndex should be 0 if the first
    node is desired.
  @param[out] nValue the integer contained in the node.
  */
int getXMLValueNoThrow(XMLNode xParent, std::string sChild,int nIndex,unsigned int &nValue);/**<
  A version of \ref getXMLValue for integers that doesn't throw an exception if no xml node is 
  found. It returns 0 if node is either empty or not found, otherwise it returns 1 to indicate
  that the node exists and that it is not empty.
  
  @param[in] xParent the parent \ref XMLNode containing the node which contains the integer value
    being seeked.
  @param[in] sChild the name of the \ref XMLNode containing the integer value being seeked.
  @param[in] nIndex the index of the child node i.e. if there is more than one node of the same
    name and the second node is desired then nIndex should be 1 and nIndex should be 0 if the first
    node is desired.
  @param[out] nValue the integer contained in the node.
  */
int getXMLValueNoThrow(XMLNode xParent, std::string sChild,int nIndex, double &dValue);/**<
  A version of \ref getXMLValue for doubles that doesn't throw an exception if no xml node is 
  found. It returns 0 if node is either empty or not found, otherwise it returns 1 to indicate
  that the node exists and that it is not empty.
  
  @param[in] xParent the parent \ref XMLNode containing the node which contains the double value
    being seeked.
  @param[in] sChild the name of the \ref XMLNode containing the double value being seeked.
  @param[in] nIndex the index of the child node i.e. if there is more than one node of the same
    name and the second node is desired then nIndex should be 1 and nIndex should be 0 if the first
    node is desired.
  @param[out] dValue the double contained in the node.
  */
int getXMLValueNoThrow(XMLNode xParent, std::string sChild,int nIndex, std::string &sValue);/**<
  A version of \ref getXMLValue for strings that doesn't throw an exception if no xml node is 
  found. It returns 0 if node is either empty or not found, otherwise it returns 1 to indicate
  that the node exists and that it is not empty.
  
  @param[in] xParent the parent \ref XMLNode containing the node which contains the string value 
    being seeked.
  @param[in] sChild the name of the \ref XMLNode containing the string value being seeked.
  @param[in] nIndex the index of the child node i.e. if there is more than one node of the same
    name and the second node is desired then nIndex should be 1 and nIndex should be 0 if the first
    node is desired.
  @param[out] sValue the string contained in the node.
  */
int getXMLValueNoThrow(XMLNode xParent, std::string sChild,int nIndex, bool &bValue);/**<
  A version of \ref getXMLValue for bools that doesn't throw an exception if no xml node is 
  found. It returns 0 if node is either empty or not found, otherwise it returns 1 to indicate
  that the node exists and that it is not empty.
  
  @param[in] xParent the parent \ref XMLNode containing the node which contains the bool value 
    being seeked.
  @param[in] sChild the name of the \ref XMLNode containing the bool value being seeked.
  @param[in] nIndex the index of the child node i.e. if there is more than one node of the same
    name and the second node is desired then nIndex should be 1 and nIndex should be 0 if the first
    node is desired.
  @param[out] sValue the bool contained in the node.
  */
int getXMLAttributeNoThrow(XMLNode xParent, std::string sAttribute, int &nValue);
int getXMLAttributeNoThrow(XMLNode xParent, std::string sAttribute, std::string &sValue);/**<
  Attempts to get the attribute named by \c sAttribute in the parent node \c xParent and returns its
  value in \c sValue. If attribute is found it returns 1, otherwise it retunrs 0.
  */
int getXMLAttributeNoThrow(XMLNode xParent, std::string sAttribute, double &dValue);/**<
  Attempts to get the attribute named by \c sAttribute in the parent node \c xParent and returns its
  value in \c dValue. If attribute is found it returns 1, otherwise it retunrs 0.
  */
XMLNode getXMLNode(XMLNode xParent, std::string sChild,int nIndex);/**<
  Returns the \c nIdex th XMLNode in \c xParent with name \c sChild. If no node is found it throws
  an exception saying that no node was found when one was expected.
  */
XMLNode getXMLNodeNoThrow(XMLNode xParent, std::string sChild,int nIndex);/**<
  returns the \c nIdex th XMLNode in \c xParent with name \c sChild. If no node is found it simply
  returns an empty node. One can test if the node is empty using XMLNode::isEmpty().
  */
XMLNode openXMLFile(std::string sFileName, std::string sStartNode);
XMLNode openXMLFileNoThrow(std::string sFileName, std::string sStartNode);
#endif
