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
    
def twoKatan(pR, KE, plain, rounds=254, faultRound=237):
  '''
  calculating the symbolic difference from two katan implementations
  '''
  #generating the systems
  origOutput, origKatan = oneKatan32enc(rounds, pR, plain)
  faultOutput, faultKatan, fault_pos, faultSys = faultKatan32(rounds, faultRound, pR, plain)
  assert len(origKatan) == len(faultKatan) and len(origOutput)==len(faultOutput), 'the equation system have not the same length'
  #for f in faultSys:
  #  print f
  fullSys=faultSys + origKatan + origOutput + faultOutput + faultKatan
  
  print 'length of the system: ', len(origKatan)

  #adding up the instances
  delta=set()
  deltaOuput=[]
  for i in xrange(len(origKatan)):
    delta.add(origKatan[i]+faultKatan[i])
  for i in xrange(32):
    deltaOuput.append(origOutput[i]+faultOutput[i])
  fullSys+=delta
  fullSys+=deltaOuput

  #non-symbolic systems (adding output to the system)
  origC = Katan32enc(plain, rounds, KE)
  faultC = Katan32fault(plain, rounds, faultRound, fault_pos, KE)

  fullSys += [ 'zD_%02d+%s'%(i, str(origC[i])) for i in range(32) ] 
  fullSys += [ 'zF_%02d+%s'%(i, str(faultC[i])) for i in range(32) ]
  return fullSys

def twoKatanVec(pR, KE, plain, rounds=254, faultRound=237, fault_vec=None, fakeVec=None):
  '''
  calculating the symbolic difference from two katan implementations
  '''
  #y = sorted([ v for v in pR.gens() if str(v)[0] == "y"], key=lambda k: str(k))
  zD = sorted([ v for v in pR.gens() if str(v)[:2] == "zD"], key=lambda k: str(k))
  zF = sorted([ v for v in pR.gens() if str(v)[:2] == "zF"], key=lambda k: str(k))
  
  if not fakeVec: fakeVec=fault_vec

  #generating the systems
  origOutput, origKatan = oneKatan32enc(rounds, pR, plain)
  faultOutput, faultKatan, faultSys = faultKatan32Vec(rounds, faultRound, pR, plain, fakeVec)
  assert len(origKatan) == len(faultKatan) and len(origOutput)==len(faultOutput), 'the equation system have not the same length'
  
  fullSys=faultSys + origKatan + origOutput + faultOutput + faultKatan
  print 'length of the system: ', len(origKatan)

  #adding up the instances
  delta=set()
  deltaOuput=[]
  for i in xrange(len(origKatan)):  
    delta.add(origKatan[i]+faultKatan[i])
  for i in xrange(32):
    deltaOuput.append(origOutput[i]+faultOutput[i])
    #deltaOuput.append(zD[i] + zF[i] + y[i])
  fullSys+=delta
  fullSys+=deltaOuput
  
  #non-symbolic systems (adding output to the system)
  origC = Katan32enc(plain, rounds, KE)
  faultC = Katan32faultVec(plain, fault_vec, rounds, faultRound, KE)

  for i in xrange(32):
    print origC[i] + faultC[i],
  print
  fullSys += [ 'zD_%02d+%s'%(i, str(origC[i])) for i in range(32) ] 
  fullSys += [ 'zF_%02d+%s'%(i, str(faultC[i])) for i in range(32) ]
  return fullSys

def oneKatan32enc(rounds, pR, plain):
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
    return_bit = key_lfsr[0]
    
    #calculating the new bit
    new_bit = key_lfsr[0]+key_lfsr[19]+key_lfsr[30]+key_lfsr[67]

    #shifting the lfsr
    key_lfsr[:79] = key_lfsr[1:]
    key_lfsr[79] = new_bit

    #return k[curRound]
    return return_bit

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
  k = sorted([ v for v in allVars if str(v)[0] == "k"], key=lambda l: str(l))
  z = sorted([ v for v in allVars if str(v)[:2] == "zD"], key=lambda k: str(k))
  global a,b
  a = sorted([ v for v in allVars if str(v)[:2] == "aD"], key=lambda k: int(str(k)[3:]))
  b = sorted([ v for v in allVars if str(v)[:2] == "bD"], key=lambda k: int(str(k)[3:]))
  
  #setting up the array for the system of equations
  global F
  F = []

  #loading the plaintext in the registers of size 13 and 19 respectively
  global L1,L2
  L1 = plain[19:]
  L2 = plain[:19]

  #LFSR for the key schedule; lsb first
  global key_lfsr
  key_lfsr = x
  
  for curRound in xrange(rounds):
    #clocking the key_lfsr twice
    key = {}
    key['a'] = key_schedule(2*curRound)
    key['b'] = key_schedule(2*curRound+1)
    update_state(key, curRound)
  
  #output generation
  output = L2+L1
  for i in range(32):
    output[i]+=z[i]
  return output, F

def faultKatan32(rounds, faultRound, pR, plain):
  """
  a research only implementation of katan32.
  INPUT:
    faultRound - the round in which the fault is inserted. No fault .. -1 (default: -1)
    plan - the plaintext in \F_2. least signicant first (default: None)
    rounds - the number of rounds
    full_key - key for the cipher in \F_2. None for symbolic (default: None)
  """
  def key_schedule(curRound):
    #producing new output
    #F.append(k[curRound] + key_lfsr[0])
    #key_lfsr[0] = k[curRound]
    return_bit = key_lfsr[0]

    #calculating the new bit
    new_bit = key_lfsr[0]+key_lfsr[19]+key_lfsr[30]+key_lfsr[67]

    #shifting the lfsr
    key_lfsr[:79] = key_lfsr[1:]
    key_lfsr[79] = new_bit

    #return k[curRound]
    return return_bit

  def update_state(key, curRound, faultRound):#updates the symbolic state
    #setting new variables
    if curRound < faultRound:
      F.append(a[curRound] + L2[0])
      F.append(b[curRound] + L1[0])
      L2[0] = a[curRound]
      L1[0] = b[curRound]
    else:
      F.append(aF[curRound] + L2[0])
      F.append(bF[curRound] + L1[0])
      L2[0] = aF[curRound]
      L1[0] = bF[curRound]

  
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
  x = sorted([ v for v in allVars if str(v)[0] == "x"], key=lambda k: str(k))
  k = sorted([ v for v in allVars if str(v)[0] == "k"], key=lambda l: str(l))
  z = sorted([ v for v in allVars if str(v)[:2] == "zF"], key=lambda k: str(k))
  global a,b,aF,bF
  a = sorted([ v for v in allVars if str(v)[:2] == "aD"], key=lambda k: int(str(k)[3:]))
  aF = sorted([ v for v in allVars if str(v)[:2] == "aF"], key=lambda k: int(str(k)[3:]))
  b = sorted([ v for v in allVars if str(v)[:2] == "bD"], key=lambda k: int(str(k)[3:]))
  bF = sorted([ v for v in allVars if str(v)[:2] == "bF"], key=lambda k: int(str(k)[3:]))
  #setting up the array for the system of equations
  global F
  F = []
  faultF = []

  #loading the plaintext in the registers of size 13 and 19 respectively
  global L1,L2
  L1 = plain[19:]
  L2 = plain[:19]

  #LFSR for the key schedule; lsb first
  global key_lfsr
  key_lfsr = x
    
  for curRound in xrange(rounds):
    if faultRound==curRound:
      #inserting fault
      fault_pos = ZZ.random_element(32)
      print 'def fp', fault_pos
      #inserting fault
      if fault_pos < 13:
        L1[fault_pos] += 1
      else:
        L2[fault_pos-13] += 1#clocking the key_lfsr twice
      
    #clocking the key_lfsr twice
    key = {}
    key['a'] = key_schedule(2*curRound)
    key['b'] = key_schedule(2*curRound+1)
    update_state(key, curRound, faultRound)
  
  #output generation
  output = L2 + L1
  for i in range(32):
    output[i]+=z[i]
  return output, F, fault_pos, faultF

def faultKatan32Vec(rounds, faultRound, pR, plain, fault_vec):
  """
  a research only implementation of katan32.
  INPUT:
    faultRound - the round in which the fault is inserted. No fault .. -1 (default: -1)
    plan - the plaintext in \F_2. least signicant first (default: None)
    rounds - the number of rounds
    full_key - key for the cipher in \F_2. None for symbolic (default: None)
  """
  def key_schedule(curRound):
    #producing new output
    #F.append(k[curRound] + key_lfsr[0])
    #key_lfsr[0] = k[curRound]
    return_bit = key_lfsr[0]

    #calculating the new bit
    new_bit = key_lfsr[0]+key_lfsr[19]+key_lfsr[30]+key_lfsr[67]

    #shifting the lfsr
    key_lfsr[:79] = key_lfsr[1:]
    key_lfsr[79] = new_bit
    
    return return_bit
    #return k[curRound]

  def update_state(key, curRound, faultRound):#updates the symbolic state
    #setting new variables
    if curRound < faultRound:
      F.append(a[curRound] + L2[0])
      F.append(b[curRound] + L1[0])
      L2[0] = a[curRound]
      L1[0] = b[curRound]
    else:
      F.append(aF[curRound] + L2[0])
      F.append(bF[curRound] + L1[0])
      L2[0] = aF[curRound]
      L1[0] = bF[curRound]

  
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
  x = sorted([ v for v in allVars if str(v)[0] == "x"], key=lambda k: str(k))
  k = sorted([ v for v in allVars if str(v)[0] == "k"], key=lambda l: str(l))
  z = sorted([ v for v in allVars if str(v)[:2] == "zF"], key=lambda k: str(k))
  global a,b,aF,bF
  a = sorted([ v for v in allVars if str(v)[:2] == "aD"], key=lambda k: int(str(k)[3:]))
  aF = sorted([ v for v in allVars if str(v)[:2] == "aF"], key=lambda k: int(str(k)[3:]))
  b = sorted([ v for v in allVars if str(v)[:2] == "bD"], key=lambda k: int(str(k)[3:]))
  bF = sorted([ v for v in allVars if str(v)[:2] == "bF"], key=lambda k: int(str(k)[3:]))
  #setting up the array for the system of equations
  global F
  F = []
  faultF = []

  #loading the plaintext in the registers of size 13 and 19 respectively
  global L1,L2
  L1 = plain[19:]
  L2 = plain[:19]

  #LFSR for the key schedule; lsb first
  global key_lfsr
  key_lfsr = x
    
  for curRound in xrange(rounds):
    if faultRound==curRound:
      #inserting fault
      #fault_pos = ZZ.random_element(32)
      #print 'def fp', fault_pos
      
      #inserting fault
      for fault_pos in fault_vec:
        if fault_pos < 13:
          L1[fault_pos] += 1
        else:
          L2[fault_pos-13] += 1
      
    #clocking the key_lfsr twice
    key = {}
    key['a'] = key_schedule(2*curRound)
    key['b'] = key_schedule(2*curRound+1)
    update_state(key, curRound, faultRound)
  
  #output generation
  output = L2 + L1
  for i in range(32):
    output[i]+=z[i]
  return output, F, faultF
