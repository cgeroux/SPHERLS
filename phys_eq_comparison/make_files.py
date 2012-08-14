#!/usr/bin/env python

import commands

#open file
f=open("functions.txt")
for line in f:
  cmd="touch "+line[0:len(line)-1]+".cpp"
  print cmd
  commands.getoutput(cmd)