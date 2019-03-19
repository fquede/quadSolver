# eqSolver
# Connects Sage with the "quadSolver" programme
# 2013-06-13
# Author: Christopher Wolf, chris/at/Christopher-Wolf.de

import sys, subprocess;

import configTrivSolver

# constants
programPath = configTrivSolver.eqSolverPath

class eqSolver(object):
  '''Sage-Wrapper to access the C++-programme "quadSolver".
     It offers methods to easily solver a (structured) quadratic equation over GF(2).
  '''
  def __init__(self):
    '''Connect the object to the C++-programme. Creats a new instance of the 
       programme for each instance of the object.
    '''
    self.process = subprocess.Popen(args=[''], executable=programPath, bufsize=int(100*1024), shell=False, \
                   stdin=subprocess.PIPE, stdout=subprocess.PIPE, close_fds=False, )
    self.nounce = randint(1,10*1000)
    self.keyTestValues = []
    self.closed = False
    
    # send the command
    if not(self.checkConnection()):
      print "Could not connect to the sub-programme"
      exit(1)

  def __del__(self):
    self.terminate()
    return


  def __sendCommand__(self, commandString, parameterString = None):
    '''Sends a command to fastVerifier. Also collects the output (1 line).
       The command string / parameter string must NOT end on '\n'
    '''
    try:
      # debug
      #print str(self.nounce) + " " + commandString + (" " + parameterString if parameterString != None else "")
      #print "write", self.nounce, commandString
      # send the command
      self.process.stdin.write(str(self.nounce) + " " + commandString)
      #if (parameterString != None): 
      #  if len(parameterString) < 25: print "  with param=", parameterString
      if parameterString: 
        self.process.stdin.write(" " + parameterString)
      self.process.stdin.write("\n"); self.process.stdin.flush()
    except: return "Error: Could not write to pipe in cubinator.__sendCommand__"
    # get the data
    fullOutput = ""; goOn = True
    output = ""
    while goOn:
      while True:
        try:
          # collect the responds
          #print "want input"
          output = self.process.stdout.readline()
          ## debug
          #print len(output), "got/%s/"%output
        except IOError: pass
        except: return "Error: Could not read from pipe in cubinator.__sendCommand__"
        if len(output) > 0: break   # O.K. - we got a line
        #sleep(1) # give readline a bit of time
      # outer while
      # debug
      #print "read", self.nounce, output
      if output == '0': return

      # extract the nounce
      strPos = 0
      while output[strPos] in ["-"] + [str(i) for i in range(10)]: strPos += 1
      resNum = int(output[:strPos]) # if strPos > 0 else 0
      if (resNum == -1): print "debug:", output; continue # was just a debug info

      # check if we need to go on
      assert (output[strPos] == '+') or (output[strPos] == ' '), "ill formed return string for %s"%output
      goOn = output[strPos] == "+"


      # skip nounce, plus leading/trailing blanks
      output = output[strPos+1:].strip()
      if resNum != self.nounce: 
        print "Error - wrong communication for nounce=%d and res=%d. Terminate"%(self.nounce, resNum)
        exit(1)
      if len(fullOutput) == 0: fullOutput = output
      else: fullOutput += "\n" + output


    # O.K. - everything seems to be fine
    self.nounce += randint(1,100) # get ready for the next run
    return fullOutput


  def checkConnection(self):
    '''Checks if we have a working handle to the fastVerifyer programme'''
    return self.__sendCommand__("checkConnection") == "check connection O.K."     


  def addEquations(self, inList, isString = False):
    '''Sends a list of equations to the system. If necessary, converts them to a list of strings first.

       * inList: list of equations over GF(2). Either polynomials or strings

       * isString: input is a list of strings

       * return: True iff everything is O.K.
    '''
    if isinstance(inList, str): outStr = inList
    else: outStr = "\n".join([ str(a) for a in inList ])
    # ensure that there is no zero-equation in the list
    outStr += "\n" + "=" # list terminator (communiction protocol)
    #print "want to send\n", outStr
    return self.__sendCommand__("addEqs", outStr) == "add eqs O.K."


  def addVariables(self, inList):
    '''Adds variables from inList to the current system. Does not add equations.

       * inList: list of variables to add (string)
    '''
    return self.addEquations([ "%s + %s"%(a,a) for a in inList ], isString=True)


  def stats(self, level=1):
    '''Returns the current statistics on the equation system.

       * level (default: 1): easy stats are level 1, more difficult statistics are level 2.
    '''
    workingString = self.__sendCommand__("stats", str(level))
    # split the string according to the fields
    firstList = workingString.split("\n")
    assert(firstList.pop() == "done stats")
    retDict = dict();
    for curLine in firstList:
      left, right = curLine.split(" ")
      retDict[left] = right
    return retDict


  def terminate(self):
    '''Terminates the underlying programme'''
    if self.closed: return 
    self.__sendCommand__("end") 
    self.closed = True
    return 


  def getStat(self, kind):
    '''Returns the corresponding statistic'''
    return self.__sendCommand__("stat" + kind) 


  def getVars(self, inList):
    '''Returns all variables from inList (string)'''
    outStr = self.__sendCommand__("getVars", " ".join(inList))  
    return outStr


  def getStatChains(self):
    '''Returns all variables from inList (string)'''
    outStr = self.__sendCommand__("statChains") 
    return outStr


  def getLinSys(self):
    '''Returns the linear full system in the solver'''
    return self.__sendCommand__("getLinSys") 


  def getQuadSys(self):
    '''Returns the linear full system in the solver'''
    return self.__sendCommand__("getQuadSys") 


  def echelonizeSparse(self, maxMonNum = 2^31):
    '''enforces echelonization of the quadratic system. 
    
       * maxMonNum: maximal number of monomials allowed in the system. Otherwise: break
    '''
    assert 0 <= maxMonNum <= 2^31, "maxMonNum out of scope: %d"%maxMonNum

    outStr = self.__sendCommand__("echelonizeSparse", str(maxMonNum)) 
    if outStr != "done echelonizeSparse": print "got error:", outStr
    return


  def echelonizeDense(self):
    '''enforces echelonization of the quadratic system'''
    outStr = self.__sendCommand__("echelonizeDense") 
    if outStr != "done echelonizeDense": print "got error:", outStr
    return

  def fullSparseEval(self):
    '''Fully evaluates the given system with the following restrictions: stops if
 
       * the system becomes larger than "maxMonNum" (non-zero entries)
 
       * an individual polynomial is increased by more than "outMax" monomials

       * the input for one fully quadratic multiplication is larger than "inMax"
    '''
    outStr = self.__sendCommand__("fullSparseEval") 
    if (outStr[:18] == "fullSEval changed "): return int(outStr[18:])
    return "Error with " + outStr

  def fullEval(self, maxMonNum=2^31, outMax=2^31, inMax=10000):
    '''Fully evaluates the given system with the following restrictions: stops if
 
       * the system becomes larger than "maxMonNum" (non-zero entries)
 
       * an individual polynomial is increased by more than "outMax" monomials

       * the input for one fully quadratic multiplication is larger than "inMax"
    '''
    assert 0 <= maxMonNum <= 2^31, "maxMonNum out of scope: %d"%maxMonNum
    assert 0 <= outMax <= 2^31, "outMax out of scope: %d"%outMax
    assert 0 <= inMax <= 2^31, "inMax out of scope: %d"%inMax
    outStr = self.__sendCommand__("fullEval", "%d %d %d"%(maxMonNum, outMax, inMax)) 
    if (outStr[:17] == "fullEval changed "): return int(outStr[17:])
    return "Error with " + outStr


  def softReduce(self, maxMonNum=2^31, maxPolyExpand=0):
    '''Reduces the equations. stops if
 
       * the system becomes larger than "maxMonNum"
 
       * an individual polynomial is increased by more than "maxPolyExpand" monomials
    '''
    assert 0 <= maxMonNum <= 2^31, "maxMonNum out of scope: %d"%maxMonNum
    assert 0 <= maxPolyExpand <= 2^31, "maxPolyExpand out of scope: %d"%maxPolyExpand
    outStr = self.__sendCommand__("softReduce", "%d %d"%(maxMonNum, maxPolyExpand)) 
    if (outStr != "done softReduce"): print "got error:", outStr
    return # everything is fine :-)


  def sl(self):
    '''enforces Sparse Linearization (SL)-step for the quadratic system'''
    outStr = self.__sendCommand__("sl") 
    if (outStr[:11] == "sl changed "): return int(outStr[11:])
    return "Error with " + outStr

  def hsl(self, density = 0.001, densityMax = 0.0016):
    '''enforces Sparse Linearization (SL)-step for the quadratic system. Creates new monomials.

       * density: density of the system

       * densityMax: if density >= densityMax, use dense echelonization

    '''
    outStr = self.__sendCommand__("hsl", "%d %d"%(density,densityMax)) 
    if (outStr[:12] == "hsl changed "): return int(outStr[12:])
    return "Error with " + outStr


def showAll():
  '''Show some functionality of the equation solver'''
  pR = BooleanPolynomialRing(160, ["k%02d"%i for i in range(80) ] + ["i%02d"%i for i in range(80) ])
  pR.inject_variables(verbose=False)

  connectionObject = eqSolver()

  print "test connection", connectionObject.checkConnection()

  print "add vars", connectionObject.addVariables(["a", "b", "c"])
  
  print connectionObject.stats();

  print "add equations", connectionObject.addEquations([" a + b", "b + c + 1", "c  + 1"], isString=True)

  print connectionObject.getVars(["a", "b", "c"]);  

  print connectionObject.stats();

  print connectionObject.getLinSys();

  return 


def testLin():
  '''Tests the linear equation solver (random sparse equations)'''

  # little initial test
  #quadSolver = eqSolver()
  ##quadSolver.addEquations(["a+b+1", "a + b", "b +c +1 ", "a+b+c", "a+d+1"])
  #quadSolver.addEquations(["1+k32+k33+0+0+0*0+0+y0001_31", "k31+k32+0+0+0*0+0+y0001_32", "k30+k31+0+0+0*0+0+y0001_33", 
  #                         "k29+k30+0+0+0*0+0+y0001_34", "1+k28+k29+0+0+0*0+0+y0001_35", "k27+k28+0+0+0*0+0+y0001_36"])
  #print quadSolver.getQuadSys()
  #print quadSolver.getLinSys()
  #print quadSolver.stats() 
  #return

  pR = PolynomialRing(GF(2), ["x%03d"%i for i in range(100) ] )
  allVars = pR.gens()
  for repCnt in range(10):
    print "try", repCnt
    curSol = dict()
    for curVar in allVars: curSol[curVar] = randint(0,1)
    #print "sol", [ "%s:%s"%(k, curSol[k]) for k in sorted(curSol.keys())]
    
    # establish a solver object
    quadSolver = eqSolver()
    quadSolver.addVariables(sorted([str(a) for a in allVars]))
    while(int(quadSolver.stats()["linEqs"]) != len(allVars)):
      # add a random equation
      curPoly = add([ allVars[randint(0,len(allVars)-1)] for terms in range(randint(3,8)) ])
      curPoly += curPoly.substitute(curSol) # fix the constant term
      assert(curPoly.substitute(curSol) == 0) # is the equation valid? (just a precaution)
      quadSolver.addEquations([str(curPoly)])
      #print
      #print "new:", curPoly
      #print quadSolver.getLinSys()
      #print quadSolver.stats()
      
    # check if the solution is correct
    fullSols = set(quadSolver.getVars([str(a) for a in allVars]).split(","))
    for curVar in allVars:
      if not("%s=%d"%(curVar, curSol[curVar])): print "Error at", curVar
    print "checked"
    print


def testQuad():
  # little initial test
  quadSolver = eqSolver()
  quadSolver.addEquations(["a+b+c", "a*b + c*d + 1", "a*b + c*a + a", "c*a + c*b +b", "c*d + a*c + a*b + c", "c*a + a*b + a", "a*c + a + b + c"])
  print quadSolver.stats(level=2)
  print "eval one:   ", quadSolver.fullEval(inMax=1)
  print quadSolver.stats(level=2)
  print "eval two:   ", quadSolver.fullEval(inMax=1)
  print "eval three: ", quadSolver.fullEval()
  print quadSolver.getLinSys()
  print quadSolver.getQuadSys()
  print quadSolver.stats(level=2)


def testQuadLarge():
  '''Tests the quadratic equation solver (random sparse equations)'''
  # little initial test
  quadSolver = eqSolver()
  inEqs = ["a*b + c*d + 1", "a*b + c*a + a", "c*a + c*b +b", 
           "c*d + a*c + a*b + c", "c*a + a*b + a", 
           "a*c + a + b + c",
           "a+b+c", "a+b+1"
  ]
  quadSolver.addEquations(inEqs)
  # ensure it is solvable
  #R = BooleanPolynomialRing(4, ["a", "b", "c", "d"], order='degneglex')
  #I = ideal([ R(a) for a in inEqs ])
  #B = I.groebner_basis()
  #print B
  #return

  #print quadSolver.getQuadSys()
  #print quadSolver.getLinSys()
  #print "eval one:   ", quadSolver.fullEval()
  #print "eval two:   ", quadSolver.fullEval()
  #print "eval three: ", quadSolver.fullEval()
  #print quadSolver.getLinSys()
  #print quadSolver.getQuadSys()
  #print quadSolver.stats()
  
  #return

  pR = PolynomialRing(GF(2), ["x%03d"%i for i in range(20) ] )
  allVars = pR.gens()
  for repCnt in range(10):
    print "try", repCnt
    sys.stdout.flush()
    curSol = dict()
    for curVar in allVars: curSol[curVar] = randint(0,1)
    #print "sol", [ "%s:%s"%(k, curSol[k]) for k in sorted(curSol.keys())]
    
    # establish a solver object
    quadSolver = eqSolver()
    quadSolver.addVariables(sorted([str(a) for a in allVars]))  # make sure all variables get the right number
    fieldEq = ideal([ a^2 + a for a in allVars ]);
    while(int(quadSolver.stats()["linRows"]) != len(allVars)):
      # add a random equation
      #curPoly = add([ allVars[randint(0,len(allVars)-1)]*allVars[randint(0,len(allVars)-1)] for terms in range(randint(3,10)) ])
      # enforce a structure: only variables in range +/-1 are allowed now
      startVar = randint(0,len(allVars)-1)
      varPos = set([startVar, (startVar+1) % len(allVars), (startVar-1) % len(allVars), randint(0,len(allVars)-1)])
      curPoly = add([randint(0,1)*allVars[a]*allVars[b] for a in varPos for b in varPos])
      # x^2 -> x
      curPoly = fieldEq.reduce(curPoly)
      curPoly += add([ allVars[randint(0,len(allVars)-1)] for terms in range(randint(0,3)) ]) # some linear terms
      curPoly += curPoly.substitute(curSol) # fix the constant term
      assert(curPoly.substitute(curSol) == 0) # is the equation valid? (just a precaution)
      quadSolver.addEquations([str(curPoly)])
      #print
      #print "new:", curPoly
      #print quadSolver.getLinSys()
      #print quadSolver.getQuadSys()
      #print quadSolver.stats()
      if randint(0,len(allVars)) == 5: 
        quadSolver.echelonizeSparse()
        print "eval"
        quadSolver.fullEval() # enforce full evaluation from time to time
        print quadSolver.stats()
        print quadSolver.getStatChains()
        sys.stdout.flush()
       
        # get the short equations
        fullEqSys = quadSolver.getQuadSys()
        for curLine in fullEqSys.split("\n"):
          if curLine.count("*") < 3: print curLine

    # some last eval
    reEval = 0
    while(int(quadSolver.stats()["quadRows"]) != 0):
      reEval += 1
      quadSolver.fullEval()
      #print quadSolver.getLinSys()
      #print quadSolver.getQuadSys()
      #print

      
    # check if the solution is correct
    fullSols = set(quadSolver.getVars([str(a) for a in allVars]).split(","))
    for curVar in allVars:
      if not("%s=%d"%(curVar, curSol[curVar])): print "Error at", curVar
    print quadSolver.stats()
    print "checked, re-eval", reEval
    print

'''def testQuadReplace():
  # little initial test
  quadSolver = eqSolver()
  # solution: a=1, b=0, c=1, d=1, e=0
  #quadSolver.addEquations(["a*b + b*c + d*e + a + b +1", "a*c + c*d", "a*b + m1", "b*c + m2", "d*e + m3"])
  #quadSolver.addEquations(["a*b + d*e + e + d + 1", "b*c + a*b + a +b + c + e", "a*c + m4", "c*d + m5"])
  #quadSolver.addEquations(["a*b + b*c", "b*c + d*e", "a+b+1","a+c+d+1", "c*d+1", "e*c + c + 1"])

  quadSolver.addEquations(["a*b + a*c + a*d + b", "b*c + b*e + e","d*e + b", "a*b + a*e + c + d", "c*d + a", "d*e + b"])
  #quadSolver.addEquations(["a*b + m1", "a*c + m2", "a*d + m3", "b*c + m4", "b*e + m5", "d*e + m6", "a*e + m7", "c*d + m8"])
  #quadSolver.addEquations(["a*b + a*m1", "a*c + a*m2", "a*d + a*m3","b*c + b*m4", "b*e + b*m5", "d*e + d*m6", "a*e + a*m7", "c*d + c*m8"])
  #quadSolver.addEquations(["a*b + b*m1", "a*c + c*m2", "a*d + d*m3","b*c + c*m4", "b*e + e*m5", "d*e + e*m6", "a*e + e*m7", "c*d + d*m8"])

  quadSolver.addEquations(["a*b", "a*c +1"])
  
  print "--input--"
  print quadSolver.getLinSys()
  print quadSolver.getQuadSys()
  
  quadSolver.echelonizeDense()
  quadSolver.fullEval()
  quadSolver.sl()
  quadSolver.echelonizeDense()
  quadSolver.fullEval()
  quadSolver.sl()
  quadSolver.echelonizeDense()
  quadSolver.fullEval()
    
  print "--output--"
  print quadSolver.getLinSys()
  print quadSolver.getQuadSys()

  return
# uncomment to test functionality of the class
testQuadReplace()
'''

