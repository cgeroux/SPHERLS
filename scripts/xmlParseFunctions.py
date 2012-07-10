import xml.etree.ElementTree as xml
import warnings
import mywarnings

def getOneChildElementText(parentElement,childElementName):
  """Returns text from the first element of name childElementName found in parentElement 
  printing a warning message if more than one found."""
  
  childElements=parentElement.findall(childElementName)
  if len(childElements)>1:
    warnings.warn("found "+str(len(childElements))+" \""+childElementName+"\" elements in parent node \""+
      parentElement.tag+"\", all but first element will be ignored.")
  elif len(childElements)==0:
    return None
  return childElements[0].text
def getOneChildElementTextToBool(parentElement,childElementName):
  """Returns true/false from the first element of name childElementName found in parentElement 
  printing a warning message if more than one found."""
  
  childElements=parentElement.findall(childElementName)
  if len(childElements)>1:
    warnings.warn("found "+str(len(childElements))+" \""+childElementName+"\" elements in parent element\""+
      parentElement.tag+"\", all but first element will be ignored.")
  elif len(childElements)==0:
    return None
  return childElements[0].text.lower() in ["true","yes"]
def getOneChildElementTextToFloat(parentElement,childElementName):
  """Returns float from the first element of name childElementName found in parentElement 
  printing a warning message if more than one found."""
  
  childElements=parentElement.findall(childElementName)
  if len(childElements)>1:
    warnings.warn("found "+str(len(childElements))+" \""+childElementName+"\" elements, all but "
      +"first element will be ignored.")
  elif len(childElements)==0:
    return None
  try:
    return float(childElements[0].text)
  except ValueError:
    raise Exception("expecting a float value in node \""+childElementName+"\" under parent node \""+
    parentElement.tag+"\" got \""+childElements[0].text+"\".")
def getOneChildElementTextToInt(parentElement,childElementName):
  """Returns an integer from the first element of name childElementName found in parentElement 
  printing a warning message if more than one found."""
  
  childElements=parentElement.findall(childElementName)
  if len(childElements)>1:
    warnings.warn("found "+str(len(childElements))+" \""+childElementName+"\" elements, all but "
      +"first element will be ignored.")
  elif len(childElements)==0:
    return None
  try:
    return int(childElements[0].text)
  except ValueError:
    raise Exception("expecting an integer value in node \""+childElementName+"\" under parent node \""+
    parentElement.tag+"\" got \""+childElements[0].text+"\".")
def getOneChildElement(parentElement,childElementName):
  """Returns the first xml element of name childElementName from parentElement printing a warning 
    saying that others are ignored if more than one found"""
  
  childElements=parentElement.findall(childElementName)
  if len(childElements)>1:
    warnings.warn("found "+str(len(childElements))+" \""+childElementName+"\" nodes, all but "
      +"first node will be ignored.")
  elif len(childElements)==0:
    return None
  return childElements[0]