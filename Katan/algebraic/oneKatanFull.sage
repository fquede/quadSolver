from polybori import *
from Katan32 import *
from utilities import *

global irregular
irregular = [1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0,
             1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0,
             1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1,
             1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0,
             0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0,
             1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0,
             0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1,
             1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1,
             1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1,
             0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1,
             0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0,
             1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1,
             0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0,
             0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1,
             0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0]
    
def oneKatan32enc(rounds, pR, numPlainBits, valueDict, toGuess=[]):
  """
  a research only implementation of katan32.
  INPUT:
    plan - the plaintext in \F_2. least signicant first (default: None)
    rounds - the number of rounds
    full_key - key for the cipher in \F_2. None for symbolic (default: None)
  """
  def key_schedule(curRound):
    #producing new output
    #F.append(k[curRound] + key_lfsr[0])
    #key_lfsr[0] = k[curRound]
    return_bit=key_lfsr[0]

    #calculating the new bit
    new_bit = key_lfsr[0]+key_lfsr[19]+key_lfsr[30]+key_lfsr[67]

    #shifting the lfsr
    key_lfsr[:79] = key_lfsr[1:]
    key_lfsr[79] = new_bit

    return return_bit
  
  def update_state_noVars(key, curRound):#updates the symbolic state
    #updating; the key is a dictionary a:k_a, b:k_b
    f_a = L1[12] + L1[7] + L1[8]*L1[5] + L1[3]*irregular[curRound] + key['a']
    f_b = L2[18] + L2[7] + L2[12]*L2[10] + L2[8]*L2[3] + key['b']
    
    # shift all registers by 1  
    L1[1:] = L1[:12]
    L1[0] = f_b
    L2[1:] = L2[:18]
    L2[0] = f_a

  def update_state(key, curRound):#updates the symbolic state
    #setting new variables
    F.append(a[curRound] + L2[0])
    F.append(b[curRound] + L1[0])
    L2[0] = a[curRound]
    L1[0] = b[curRound]

    #updating; the key is a dictionary a:k_a, b:k_b
    f_a = L1[12] + L1[7] + L1[8]*L1[5] + L1[3]*irregular[curRound] + key['a']
    f_b = L2[18] + L2[7] + L2[12]*L2[10] + L2[8]*L2[3] + key['b']
    
    # shift all registers by 1  
    L1[1:] = L1[:12]
    L1[0] = f_b
    L2[1:] = L2[:18]
    L2[0] = f_a

  #loading the variables
  allVars=pR.gens()
  p = sorted([ v for v in allVars if str(v)[0] == "p"], key=lambda k: str(k))
  x = sorted([ v for v in allVars if str(v)[0] == "x"], key=lambda k: str(k))
  if toGuess!=[]:
    for i in toGuess:
      x[i] = GF(2)(valueDict[pR('x%02d'%i)])
  #k = sorted([ v for v in allVars if str(v)[0] == "k"], key=lambda l: str(l))
  c = sorted([ v for v in allVars if str(v)[0] == "c"], key=lambda k: str(k))
  global a,b
  a = sorted([ v for v in allVars if str(v)[0] == "a"], key=lambda k: int(str(k)[3:]))
  b = sorted([ v for v in allVars if str(v)[0] == "b"], key=lambda k: int(str(k)[3:]))
  
  #setting up the array for the system of equations
  global F
  F = []

  #loading the plaintext in the registers of size 13 and 19 respectively
  global L1,L2
  L1 = [0 for _ in xrange(13)]
  L2 = [0 for _ in xrange(19)]
  for i in xrange(numPlainBits):
    plainName = "pD_%02d"%i
    if i < 13: L1[12-i] = pR(plainName)
    else: L2[31-i] = pR(plainName)

  #LFSR for the key schedule; lsb first
  global key_lfsr
  key_lfsr = x
  
  rounds_woVars = 0
  #if toGuess!=[]:
  #  rounds_woVars=17
  #  for curRound in xrange(rounds_woVars):
  #    #clocking the key_lfsr twice
  #    key = {}
  #    key['a'] = key_schedule(2*curRound)
  #    key['b'] = key_schedule(2*curRound+1)
  #    update_state_noVars(key, curRound)
  #  
  #  output = L2 + L1
  #  #for i in range(32):
  #  #  print output[i]
  #  for i in range(12):
  #    F.append(b[rounds_woVars-i-1] + L1[i+1])
  #    L1[i+1] = b[rounds_woVars-i-1]
  #  for i in range(18):
  #    F.append(a[rounds_woVars-i-1] + L2[i+1])
  #    L2[i+1] = a[rounds_woVars-i-1]
  
  for curRound in xrange(rounds_woVars, rounds):
    #clocking the key_lfsr twice
    key = {}
    key['a'] = key_schedule(2*curRound)
    key['b'] = key_schedule(2*curRound+1)
    update_state(key, curRound)
  
  #output generation
  output = L2 + L1
  for i in range(32):
    output[i]+=c[i]
  return F+output

def oneKatan32dec(rounds, pR, numPlainBits, valueDict, toGuess=[]):
  """
  a research only implementation of katan32.
  INPUT:
    plan - the plaintext in \F_2. least signicant first (default: None)
    rounds - the number of rounds
    full_key - key for the cipher in \F_2. None for symbolic (default: None)
  """
  def update_state(key, curRound, counter):#updates the symbolic state
    #setting new variables
    f_a = L2[0]
    f_b = L1[0]

    # shift all registers by 1  
    L1[:12] = L1[1:]
    L2[:18] = L2[1:]
    
    if curRound >= 18:
      L1[12] = b[counter]
      L2[18] = a[counter]
      F.append(b[counter] + f_a + L1[7] + L1[8]*L1[5] + L1[3]*irregular[curRound] + key[2*curRound])
      F.append(a[counter] + f_b + L2[7] + L2[12]*L2[10] + L2[8]*L2[3] + key[2*curRound+1])
    elif curRound >= 12 and curRound < 18:
      L1[12] = b[counter]
      L2[18] = h[curRound]
      F.append(b[counter] + f_a + L1[7] + L1[8]*L1[5] + L1[3]*irregular[curRound] + key[2*curRound])
      F.append(h[curRound] + f_b + L2[7] + L2[12]*L2[10] + L2[8]*L2[3] + key[2*curRound+1])
    elif curRound < 12:
      L1[12] = g[curRound]
      L2[18] = h[curRound]
      F.append(g[curRound] + f_a + L1[7] + L1[8]*L1[5] + L1[3]*irregular[curRound] + key[2*curRound])
      F.append(h[curRound] + f_b + L2[7] + L2[12]*L2[10] + L2[8]*L2[3] + key[2*curRound+1])

  #loading the variables
  allVars=pR.gens()
  p = sorted([ v for v in allVars if str(v)[0] == "p"], key=lambda k: str(k))
  c = sorted([ v for v in allVars if str(v)[0] == "c"], key=lambda k: str(k))
  x = sorted([ v for v in allVars if str(v)[0] == "x"], key=lambda k: str(k))
  if toGuess!=[]:
    for i in toGuess:
      x[i] = GF(2)(valueDict[pR('x%02d'%i)])
  #k = sorted([ v for v in allVars if str(v)[0] == "k"], key=lambda l: str(l))
  global a,b
  a = sorted([ v for v in allVars if str(v)[0] == "a"], key=lambda k: int(str(k)[3:]), reverse=True)
  a=a[18:]
  b = sorted([ v for v in allVars if str(v)[0] == "b"], key=lambda k: int(str(k)[3:]), reverse=True)
  b=b[12:]
  g = sorted([ v for v in allVars if str(v)[0] == "g"], key=lambda k: int(str(k)[3:]), reverse=True)
  h = sorted([ v for v in allVars if str(v)[0] == "h"], key=lambda k: int(str(k)[3:]), reverse=True)
  
  #setting up the array for the system of equations
  global F
  F = []

  #loading the plaintext in the registers of size 13 and 19 respectively
  global L1,L2
  L1 = c[19:]
  L2 = c[:19]
  
  #LFSR for the key schedule; lsb first
  global key
  key = [0 for _ in xrange(2*rounds + 80)]
  key[:80] = x
  for i in xrange(2*rounds):
    #F.append(k[i] + key[i])
    #key[i] = k[i]
    key[80+i] = key[i] + key[i+19] + key[i+30] + key[i+67]

  #rounds backwards
  counter = 0
  for curRound in xrange(rounds-1,-1,-1):
    update_state(key, curRound, counter)
    counter += 1

  #output generation
  inputP = L2+L1
  for i in xrange(numPlainBits):
    plainName = "pD_%02d"%i
    inputP[31-i] += pR(plainName)
  
  return F+inputP

#################################################################################
def specialInstance(plainCounter, eqSys, outSys, rounds, varValues, plainVarNum, noOutput):
  polyRing = parent(varValues.keys()[0])
  # transfer the plain text counter to a setting of the iv variables
  plainDict = {}
  workingVal = plainCounter
  for i in range(plainVarNum):
    plainDict['pD_%02d'%i] = str(1) if ((workingVal % 2) == 1) else str(0)
    varValues[polyRing('pD_%02d'%i)] = 1 if ((workingVal % 2) == 1) else 0
    workingVal //= 2

  print "plainDict: ",prettyDictStr(plainDict)
  
  workingSys = deepcopy(eqSys)
  workingOutSys = deepcopy(outSys)

  # generate the cipher text
  if not(noOutput):
    #creating the plain text
    plain = [0 for _ in xrange(32)]
    for i in xrange(plainVarNum):
      plain[31-i] = GF(2)(varValues[polyRing('pD_%02d'%i)])

    #key setup
    key = [0 for _ in xrange(80)]
    for i in xrange(80):
      key[i] = GF(2)(varValues[polyRing('x%02d'%i)])
  
    cipher = Katan32enc(plain,rounds,key)
    cipherDict = {}
    for i in xrange(32):
      cipherDict['cD_%02d'%i] = str(cipher[i])

    for curCipher in cipherDict.keys(): # slow but working
      workingSys = [ a.replace(curCipher, cipherDict[curCipher]) for a in workingSys ]
      workingOutSys = [ a.replace(curCipher, cipherDict[curCipher]) for a in workingOutSys ]
  
  # replace all plain variables by their corresponding value
  # replace all plain&output variables by their corresponding value
  
  #workingSys=Sequence(workingSys); workingSys = workingSys.subs(plainDict); workingSys = workingSys.subs(cipherDict)
  #workingOutSys=Sequence(workingOutSys); workingOutSys = workingOutSys.subs(plainDict); workingOutSys = workingOutSys.subs(cipherDict)
  
  for curPlain in plainDict.keys(): # slow but working
    workingSys = [ a.replace(curPlain, plainDict[curPlain]) for a in workingSys ]
    workingOutSys = [ a.replace(curPlain, plainDict[curPlain]) for a in workingOutSys ]
  
  # move the system of equations to an array of strings
  workingSys = [ str(a).replace(" ", "") for a in workingSys ]
  workingOutSys = [ str(a).replace(" ", "") for a in workingOutSys ]

  #adding output to the system
  workingSys+=workingOutSys

  # replace all "D" by ivCounter (4 digits)
  workingSys = [ s.replace("D", '%04d'%plainCounter) for s in workingSys]

  return workingSys

#parallel data collector
@parallel
def katanOutVals(instance, infoList):
  print "instance", instance,
  sys.stdout.flush()
  plainVarNum, varDict, rounds = infoList
  outVals = []
  # set plain variables
  instanceNum = int(instance)
  for i in range(plainVarNum): 
    varDict["p%02d"%i] = 1 if ((instanceNum % 2) == 1) else 0
    instanceNum //= 2

  #creating the plain text
  plain = [0 for _ in xrange(32)]
  for i in xrange(plainVarNum):
    plain[i] = GF(2)(varDict['p%02d'%i])

  #key setup
  key = [0 for _ in xrange(80)]
  for i in xrange(80):
    key[i] = GF(2)(varDict['x%02d'%i])
  
  cipher = Katan32enc(plain,rounds,key)
  #init the state of katan
  for i in range(32):
    outVals.append('z%04d_%02d+%s'%(int(instance), i, cipher[i]))

  return outVals

#parallel data collector
@parallel
def katanData(instance, infoList):
  def key_schedule(curRound):
    #producing new output
    round_key = key_lfsr[0]

    #calculating the new bit
    new_bit = key_lfsr[0]+key_lfsr[19]+key_lfsr[30]+key_lfsr[67]

    #shifting the lfsr
    key_lfsr[:79] = key_lfsr[1:]
    key_lfsr[79] = new_bit

    return round_key

  def update_state(key, curRound):#updates the symbolic state
    #updating; the key is a dictionary a:k_a, b:k_b
    f_a = L1[12] + L1[7] + L1[8]*L1[5] + L1[3]*irregular[curRound] + key['a']
    f_b = L2[18] + L2[7] + L2[12]*L2[10] + L2[8]*L2[3] + key['b']
    
    # shift all registers by 1  
    L1[1:] = L1[:12]
    L1[0] = f_b
    L2[1:] = L2[:18]
    L2[0] = f_a

  print "instance", instance,
  sys.stdout.flush()
  plainVarNum, varDict, allVars = infoList
  allRounds = [ int(a.split("_")[1]) for a in allVars ]
  minRound = min(allRounds); maxRound = max(allRounds) 
  dataEq = []
  # set plain variables
  instanceNum = int(instance)
  for i in range(plainVarNum): 
    varDict["p%02d"%i] = 1 if ((instanceNum % 2) == 1) else 0
    instanceNum //= 2
  
  #init the state of katan
  #loading the plaintext in the registers of size 13 and 19 respectively
  L1 = [0 for _ in xrange(13)]
  L2 = [0 for _ in xrange(19)]
  for i in xrange(plainVarNum):
    plainName = "p%02d"%i
    if i < 13: L1[12-i] = GF(2)(varDict[plainName])
    else: L2[31-i] = GF(2)(varDict[plainName])

  #LFSR for the key schedule; lsb first
  key_lfsr = [0 for _ in xrange(80)]
  for i in xrange(80):
    key_lfsr[i] = GF(2)(varDict['x%02d'%i])

  #extracting the variables
  for curRound in xrange(minRound):
    #clocking the key_lfsr twice
    key = {}
    key['a'] = key_schedule(2*curRound)
    key['b'] = key_schedule(2*curRound+1)
    update_state(key, curRound)
  for curRound in range(minRound, maxRound+1):
    if curRound not in allRounds: 
      key = {}
      key['a'] = key_schedule(2*curRound)
      key['b'] = key_schedule(2*curRound+1)
      update_state(key, curRound)
    else:
      curVars = [ a for a in allVars if int(a.split('_')[1]) == curRound ]
      key = {}
      key['a'] = key_schedule(2*curRound)
      key['b'] = key_schedule(2*curRound+1)
      for curVar in curVars:
        if curVar.startswith('a'): dataEq.append(curVar + '+%d'%L2[0])
        if curVar.startswith('b'): dataEq.append(curVar + '+%d'%L1[0])
      f_a = L1[12] + L1[7] + L1[8]*L1[5] + L1[3]*irregular[curRound] + key['a']
      f_b = L2[18] + L2[7] + L2[12]*L2[10] + L2[8]*L2[3] + key['b']
    
      # shift all registers by 1  
      L1[1:] = L1[:12]
      L1[0] = f_b
      L2[1:] = L2[:18]
      L2[0] = f_a

  return dataEq
