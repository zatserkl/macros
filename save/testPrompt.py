#! /usr/bin/python

import subprocess
import sys
import os

def extractCommands(commandsFileName):
  commandsFile=open(commandsFileName)
  return "".join(commandsFile.readlines())

def interrogatePrompt(commandsFileName):
  """
  Open root and pipe commands to the prompt.
  Verify the output once the commands are finished.
  """

  commands = extractCommands(commandsFileName)

  r=subprocess.Popen(["root","-l","-b"],
                     shell=False,
                     stdin=subprocess.PIPE,
                     stderr=subprocess.PIPE,
                     stdout=subprocess.PIPE)

  out,err = r.communicate(commands)
  if len(out)>0:
    print "*** Stdout:\n"+out
  if len(err)>0:
    print "*** Stderr:\n"+err


if __name__ == "__main__":
  if len(sys.argv) != 2:
    print "Usage %s commandsFileName" %os.path.basename(__file__)
    sys.exit(1)
  interrogatePrompt(sys.argv[1])
