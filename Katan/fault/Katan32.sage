from polybori import *

def Katan32enc(plain, rounds=254, full_key=None):
  """
  a research only implementation of katan32.
  INPUT:
    plan - the plaintext in \F_2. least signicant first (default: None)
    rounds - the number of rounds
    full_key - key for the cipher in \F_2. None for symbolic (default: None)
  """
  def key_schedule():
    #producing new output
    round_key = key_lfsr[0]

    #calculating the new bit
    new_bit = key_lfsr[0]+key_lfsr[19]+key_lfsr[30]+key_lfsr[67]

    #shifting the lfsr
    key_lfsr[:79] = key_lfsr[1:]
    key_lfsr[79] = new_bit

    return round_key

  def update_state(key, curRound):#updates the symbolic state
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

    #updating; the key is a dictionary a:k_a, b:k_b
    f_a = L1[12] + L1[7] + L1[8]*L1[5] + L1[3]*irregular[curRound] + key['a']
    f_b = L2[18] + L2[7] + L2[12]*L2[10] + L2[8]*L2[3] + key['b']
    	
    # shift all registers by 1	
    L1[1:] = L1[:12]
    L1[0] = f_b
    L2[1:] = L2[:18]
    L2[0] = f_a

  #loading the plaintext in the registers of size 13 and 19 respectively
  #print "plaintext: ", ''.join([str(i) for i in plain])
  global L1,L2
  L1 = plain[19:]
  L2 = plain[:19]

  #setting up the BooleanPolynomialRing and the LFSR of the key schedule; lsb first
  if full_key: #non-symbolic
    x = deepcopy(full_key)
  else: #symbolic
    R = BooleanPolynomialRing(80,'x')
    x = list(R.gens())

  #print "key: ", ''.join([str(i) for i in x])
  global key_lfsr
  key_lfsr = x
  
  for curRound in xrange(rounds):
    #clocking the key_lfsr twice
    key = {}
    key['a'] = key_schedule()
    key['b'] = key_schedule()
    update_state(key, curRound)
  
  #output are the two registers
  output = L2+L1
  #print "ciphertext: ", ''.join([str(i) for i in output])
  return output

def Katan32fault(plain, rounds=254, faultRound=190, fault_pos=17, full_key=None):
  """
  a research only implementation of katan32.
  INPUT:
    plan - the plaintext in \F_2. least signicant first (default: None)
    rounds - the number of rounds
    full_key - key for the cipher in \F_2. None for symbolic (default: None)
  """
  def key_schedule():
    #producing new output
    round_key = key_lfsr[0]

    #calculating the new bit
    new_bit = key_lfsr[0]+key_lfsr[19]+key_lfsr[30]+key_lfsr[67]

    #shifting the lfsr
    key_lfsr[:79] = key_lfsr[1:]
    key_lfsr[79] = new_bit

    return round_key

  def update_state(key, curRound, faultRound):#updates the symbolic state
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
    
    #updating; the key is a dictionary a:k_a, b:k_b
    f_a = L1[12] + L1[7] + L1[8]*L1[5] + L1[3]*irregular[curRound] + key['a']
    f_b = L2[18] + L2[7] + L2[12]*L2[10] + L2[8]*L2[3] + key['b']
      
    # shift all registers by 1  
    L1[1:] = L1[:12]
    L1[0] = f_b
    L2[1:] = L2[:18]
    L2[0] = f_a

  #loading the plaintext in the registers of size 13 and 19 respectively
  #print "plaintext: ", ''.join([str(i) for i in plain])
  global L1,L2
  L1 = plain[19:]
  L2 = plain[:19]

  #setting up the BooleanPolynomialRing and the LFSR of the key schedule; lsb first
  if full_key: #non-symbolic
    x = deepcopy(full_key)
  else: #symbolic
    R = BooleanPolynomialRing(80,'x')
    x = list(R.gens())

  #print "key: ", ''.join([str(i) for i in x])
  global key_lfsr
  key_lfsr = x
  
  for curRound in xrange(rounds):
    if faultRound==curRound:
      #inserting fault
      if fault_pos < 13:
        L1[fault_pos] += 1
      else:
        L2[fault_pos-13] += 1#clocking the key_lfsr twice
    key = {}
    key['a'] = key_schedule()
    key['b'] = key_schedule()
    update_state(key, curRound, faultRound)
  
  #output are the two registers
  output = L2+L1
  #print "ciphertext: ", ''.join([str(i) for i in output])
  return output


def Katan32faultVec(plain, fault_vec, rounds=254, faultRound=190, full_key=None):
  """
  a research only implementation of katan32.
  INPUT:
    plan - the plaintext in \F_2. least signicant first (default: None)
    rounds - the number of rounds
    full_key - key for the cipher in \F_2. None for symbolic (default: None)
  """
  def key_schedule():
    #producing new output
    round_key = key_lfsr[0]

    #calculating the new bit
    new_bit = key_lfsr[0]+key_lfsr[19]+key_lfsr[30]+key_lfsr[67]

    #shifting the lfsr
    key_lfsr[:79] = key_lfsr[1:]
    key_lfsr[79] = new_bit

    return round_key

  def update_state(key, curRound, faultRound):#updates the symbolic state
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
    
    #updating; the key is a dictionary a:k_a, b:k_b
    f_a = L1[12] + L1[7] + L1[8]*L1[5] + L1[3]*irregular[curRound] + key['a']
    f_b = L2[18] + L2[7] + L2[12]*L2[10] + L2[8]*L2[3] + key['b']
      
    # shift all registers by 1  
    L1[1:] = L1[:12]
    L1[0] = f_b
    L2[1:] = L2[:18]
    L2[0] = f_a

  #loading the plaintext in the registers of size 13 and 19 respectively
  #print "plaintext: ", ''.join([str(i) for i in plain])
  global L1,L2
  L1 = plain[19:]
  L2 = plain[:19]

  #setting up the BooleanPolynomialRing and the LFSR of the key schedule; lsb first
  if full_key: #non-symbolic
    x = deepcopy(full_key)
  else: #symbolic
    R = BooleanPolynomialRing(80,'x')
    x = list(R.gens())

  #print "key: ", ''.join([str(i) for i in x])
  global key_lfsr
  key_lfsr = x
  
  for curRound in xrange(rounds):
    if faultRound==curRound:
      #inserting fault
      for fault_pos in fault_vec:
        if fault_pos < 13:
          L1[fault_pos] += 1
        else:
          L2[fault_pos-13] += 1
    #clocking the key_lfsr twice
    key = {}
    key['a'] = key_schedule()
    key['b'] = key_schedule()
    update_state(key, curRound, faultRound)
  
  #output are the two registers
  output = L2+L1
  #print "ciphertext: ", ''.join([str(i) for i in output])
  return output


def Katan32dec(cipher, rounds=254, full_key=None):
  """
  a research only implementation of katan32.
  INPUT:
    plan - the plaintext in \F_2. least signicant first (default: None)
    rounds - the number of rounds
    full_key - key for the cipher in \F_2. None for symbolic (default: None)
  """
  def update_state(curRound):#updates the symbolic state
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
    #updating; the key is a dictionary a:k_a, b:k_b
    f_a = L2[0]
    f_b = L1[0]

    # shift all registers by 1	
    L1[:12] = L1[1:]
    L2[:18] = L2[1:]

    L1[12] = f_a + L1[7] + L1[8]*L1[5] + L1[3]*irregular[curRound] + key[2*curRound]
    L2[18] = f_b + L2[7] + L2[12]*L2[10] + L2[8]*L2[3] + key[2*curRound+1]

  #loading the plaintext in the registers of size 13 and 19 respectively
  print "ciphertext: ", ''.join([str(i) for i in cipher])
  global L1,L2
  L1 = cipher[19:]
  L2 = cipher[:19]
  
  #setting up the BooleanPolynomialRing and the LFSR of the key schedule; lsb first
  if full_key: #non-symbolic
    x = deepcopy(full_key)
  else: #symbolic
    R = BooleanPolynomialRing(80,'x')
    x = R.gens()
  print "key: ", ''.join([str(i) for i in x])
  
  #expand the key
  global key
  key = [0 for _ in xrange(2*rounds)]
  key[:80] = x
  for i in xrange(80,2*rounds):
  	key[i] = key[i-80] + key[i-61] + key[i-50] + key[i-13]
  
  for curRound in xrange(rounds-1,-1,-1):
    update_state(curRound)
  
  #output are the two registers
  output = L2 + L1
  print "plaintext: ", ''.join([str(i) for i in output])
  return output
