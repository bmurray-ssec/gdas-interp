'''
Created on Apr 27, 2011

@author: nickb
'''


class File():
  """
  Designed to reside in a Variable's self.inFile and self.outFile, these "File"
  objects provide a common interface for doing our routine functions, namely:
  
  Read/write a variable's data, variable's attributes, or file's attributes.
  """
  
  def __init__(self, filename):
    self.filename = filename
  
  
  def readData(self, varName):
    pass
    
  def readAttribute(self, varName, attrName):
    pass
  
  
  def writeData(self, writeVarName, data):
    pass
  
  def writeAttribute(self, writeVarName, attrName, attrVal):
    pass
  
  def writeVariable(self, Variable, writeDir=None, writeVarName=None):
    pass