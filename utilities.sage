#some utilities optianaly used by trivSolver and other modules

################################################################
# helper functions for strings
################################################################

##########################################################################

# collect our possible characters for shortNumbers
__num2letter = []
__num2letter.extend([str(i) for i in range(10) ])
__num2letter.extend([chr(ord('a') + i) for i in range(26) ])
__num2letter.extend([chr(ord('A') + i) for i in range(26) ])


def shortNumber(inNum, minLen = 1):
  '''Translates a number into a short string'''
  assert minLen > 0, "minLin is too small with " + str(minLen)


  outStr = ""; baseNum = len(__num2letter)
  while (len(outStr) < minLen) or (inNum > 0):
    outStr = __num2letter[inNum % baseNum] + outStr
    inNum //= baseNum
  
  return outStr


def trimStr(a):
  '''Returns a string without whitespace'''
  return str(a).strip().replace(" ", "")


def getVars(inLine):
  '''Returns the variables of a given ANF equations as a set'''
  if isinstance(inLine, str): 
    assert inLine.find(" ") == -1, "got a space in >%s<"%inLine
    a = set(flatten([ m.split("*") for m in inLine.split('+') ]))
  else:   
    a = set(flatten([ m.split("*") for m in inLine ]))
  a.discard("0"); a.discard("1"); a.discard("")
  return a


def strDeg(inPoly):
  '''Returns the degree of an input polynomial'''
  if inPoly == "0": return -1
  if inPoly == "1": return 0
  return max([ a.count("*")+1 for a in inPoly.split("+") if a != "1" ])

def varname2key(varname):
  '''Converts a given variable name to a key for sort'''
  if varname[0] == 'x': return (601, "x%02d"%(79-int(varname[1:])), varname)
  if varname[0] == 'k': return (600, int(varname[1:]), varname)
  if varname[0] == 'z': return (600, varname, int(varname.split('_')[1]))
  if varname[0] == 'y': return (603, varname, int(varname[1:]))

  assert (varname[0] in ["a", "b"]), "Wrong varname >%s<"%varname
  prefix, varRound = varname.split("_")
  varRound = int(varRound)
  if prefix[1]=='D': varInstance=1000
  else:varInstance = int(prefix[1:])
  #if varInstance.isdigit(): varInstance = int(varInstance)
  
  if prefix[0] == 'a': return (varInstance+1, -1, -varRound)
  if prefix[0] == 'b': return (varInstance+1, -2, -varRound)

def highLowRound(fullEq, high = False):
  '''Determines for a given equation its highest / lowest round (except for x-variables)''' 
  assert fullEq.find(" ") < 0, "Ill formed equation >%s<"%(fullEq)
  allRounds = [ varname2key(a)[0] for a in getVars(fullEq) if a[0] != "x" ]
  if len(allRounds) == 0: return None
  return max(allRounds) if high else min(allRounds)


##########################################################################

def str2setSet(inLine):
  '''Given a string, makes the polynomial into a set of frozen sets'''
  if isinstance(inLine, str): tmpList = inLine.split("+")
  else: tmpList = inLine
  return set([ frozenset(a.split('*')) for a in tmpList])


def setSet2str(inSet):
  return "+".join([ "*".join(a) for a in inSet ])

################################################################
# helper functions for lists of variables
################################################################

##########################################################################
def extractVars(ll,strToEx):
  #extracting a list of variables out of the list of equations ll which starts with strToEx
  #extracting polys
  listVars = [ x for x in ll if strToEx in x]
  #extracting Vars with strToEx
  exVars = set()
  for curEq in list:
    allEqVars = []
    allEqVars.extend(curEq.split(" + "))
    for var in allEqVars:
      if var[0] == strToEx:
        exVars.add(var)
  return sorted(exVars, key=lambda x: (int(x.split("_")[1]), x))


################################################################
# i/o for time and memory
################################################################

#################################################################################
# time information
def prettyTimeStr(inSecs):
  '''Gets seconds as a float, outputs a string up to days'''
  seconds = float(inSecs)
  minutes = floor(seconds) // 60; seconds -= minutes *60
  hours = minutes // 60; minutes -= hours *60
  days = hours // 24; hours -= days *24

  outStr = ""
  if days != 0: outStr += str(days) + "d ";
  if hours != 0: outStr += str(hours) + "h ";
  if minutes != 0: outStr += str(minutes) + "min ";
  outStr += "%.2f"%(seconds) + "sec";

  return outStr

#################################################################################
# memory information
def prettyMemoryStr(inBytes):
  '''Gets bytes as a float, outputs a string up GB'''
  curBytes = float(inBytes)
  unit = ""
  while (curBytes > 1000):
    if (unit == " GB"): unit = " ??"
    if (unit == " MB"): unit = " GB"
    if (unit == " kB"): unit = " MB"
    if (unit == ""): unit = " kB"
    curBytes /= 1000
  return "%.2f"%(curBytes) + unit


_proc_status = '/proc/%d/status' % os.getpid()

_scale = {'kB': 1024.0, 'mB': 1024.0*1024.0,
          'KB': 1024.0, 'MB': 1024.0*1024.0}

def _VmB(VmKey):
    '''Private.
    '''
    global _proc_status, _scale
     # get pseudo file  /proc/<pid>/status
    try:
        t = open(_proc_status)
        v = t.read()
        t.close()
    except:
        return 0.0  # non-Linux?
     # get VmKey line e.g. 'VmRSS:  9999  kB\n ...'
    i = v.index(VmKey)
    v = v[i:].split(None, 3)  # whitespace
    if len(v) < 3:
        return 0.0  # invalid format?
     # convert Vm value to bytes
    return float(v[1]) * _scale[v[2]]

def getMemory(since=0.0):
    '''Return memory usage in bytes.
    '''
    return _VmB('VmSize:') - since

def getResidentMemory(since=0.0):
    '''Return resident memory usage in bytes.
    '''
    return _VmB('VmRSS:') - since

def getStacksize(since=0.0):
    '''Return stack size in bytes.
    '''
    return _VmB('VmStk:') - since

##########################################################################

def prettyDictStr(inDict):
  '''Gets a dictionary and makes it into a string (sorted by key)'''
  outStr = ''
  first = True
  for curKey in sorted(inDict.keys()):
    if not(first): outStr += ", "
    outStr += "%s: %s"%(curKey, inDict[curKey])
    first = False

  return outStr

################################################################
# for reading and writing files
################################################################
def write_file(fileName, writeString):
  #writes a file into fileName with string
  f = open(orderFileName, "w")
  print >>f, "# File created at", time.strftime("%Y-%m-%d %H:%M:%S")  
  print >>f, "# all parameters were:", prettyDictStr(vars(args))
  print >>f
  print >>f, writeString
  f.close()
  
