#painting a matrix for the system inSys and stores it in name
from utilities import *

def paintPic(paint,inSys,name):
  # painting pictures of the monomial matrix in the quadSolver
  #for paint = 1 returns a picture of the monomial structure
  #for paint = 2 returns a picture of the variable structure paint = 1
  #for paint = 3 returns a picture of the monomials and variables that are used in the system

  #building the Ring for the matrix
  inSys = inSys.replace(" ","").split("\n")
  inSys.pop()
  termSys = flatten( [ a.split('+') for a in inSys ] )
  varSys = set(flatten( [ a.split("*") for a in termSys ])) 

  #sorting the variables
  varSys = varSys.difference(['0','1'])
  varSys = sorted(varSys, key=varname2key)
  picR = BooleanPolynomialRing(len(varSys), varSys, order = 'degneglex')
  
  #getting the system to PolyBoRi
  print "\nprocessing " + str(len(inSys)) + " Elements.."
  matSys = Sequence([picR(f) for f in inSys])
  
  #paint monomial structure
  monMat, monSem = matSys.coefficient_matrix()
  print "Matrix has size %d x %d"%(monMat.nrows(),monMat.ncols()),
  print "and rank ", monMat.rank()
  matrix_plot(monMat).save('monMat_' + str(name) + '.png')
  o = open("monSem_" + str(name) + ".txt", "w")
  for i in xrange(monSem.nrows()):
    o.write("%d: "%(i) + str(monSem[i][0]) + "\n")
  o.close()
