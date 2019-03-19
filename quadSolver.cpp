// Linear and quadratic polynomials in C++
// Good balance between readability & efficiency
//
// Authors: Christopher Wolf, chris/at/Christopher-Wolf.de, Frank Quedenfeld, frank.quedenfeld/at/googlemail.com
// 2013-05-21

// no assertions: uncomment line below; assert: leave the line commented
#define NDEBUG 

// basic includes
#include <algorithm>
#include <climits>
#include <stdio.h>
#include <assert.h>
#include <string>
#include <stdlib.h>
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */

// io
#include <fstream>
#include <iostream>
#include <sstream>

// stl data types
#include <array>
#include <bitset>
#include <deque>
#include <map>
#include <queue>
#include <set>
#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// M4RI
#include <m4ri/m4ri_config.h>
#define NDEBUG
//#undef NDEBUG
#include <m4ri/m4ri.h>

using namespace std;

typedef unsigned int uint32;
typedef vector<uint32> monVecType;


// special compare operator for variables
bool cmpVars(uint32 a, uint32 b) {  
  return a > b;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
// copies the (normalized) version of src to dest. Overwrites "dest" in the process
// note: only linear monomials "xor list of variables"
void normalizeMons(monVecType &dest, monVecType &src) {
  dest.clear();  // delete dest
  if (src.size() == 0) return; // nothing left to do
  sort(src.begin(), src.end(), cmpVars);  // get src in shape
  size_t curPos = 0;
  // delete doubles in src
  while (curPos < src.size()-1) {
    if (src[curPos] == src[curPos+1]) {
      // two doubles - delete!
      curPos += 2;
    } else {
      // no doubles - we can copy an element to workingPoly
      dest.push_back(src[curPos]);
      curPos++;
    }
  }
  // make sure we do not forget the end
  if (curPos < src.size()) dest.push_back(src[curPos]);
}

// as above, but normalizes in place. Slightly less efficient, but avoids copying
void normalizeMons(monVecType &dest) {
  if (dest.size() <= 1) return; // nothing to do here
  sort(dest.begin(), dest.end(), cmpVars);  // get dest in shape
  size_t readPos = 0;
  size_t writePos = 0;
  // delete doubles in src
  while (readPos < dest.size()-1) {
    if (dest[readPos] == dest[readPos+1]) {
      // two doubles - delete!
      readPos += 2;
    } else {
      // no doubles - we can copy an element to workingPoly
      dest[writePos] = dest[readPos];
      readPos++; writePos++;
    }
  }
  // make sure we do not forget the end
  if (readPos < dest.size()) { dest[writePos] = dest[readPos]; writePos++; }
  
  // delete the last (readPos-writePos) elements
  if (readPos > writePos) dest.resize(writePos);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
const unsigned char maxPolyDeg = 2;  // how many variables at most in a monomial?

// s[1] > s[0]
// *except* s[1] == 0. In this case, s[0] contains the linear term. 
// s[1]=s[0]=0 is the one monomial ("+1")
typedef struct {
  uint32 s[maxPolyDeg]; // array inside of a struct---needed for C++
  // mimic array access to monom_t
        uint32 &operator[] (size_t idx)             { assert(idx < maxPolyDeg); return s[idx]; }
  const uint32 &operator[] (const size_t idx) const { assert(idx < maxPolyDeg); return s[idx]; }
} monom_t;


typedef vector<monom_t> monVecQuad_t;


// checks if a given monomial is sorted
bool isSortedMonom(const monom_t &inMon) {
  for (auto i = 1; i < maxPolyDeg; i++) {
    if (inMon[i-1] == 0) continue; // ignore leading zeros
    if (inMon[i-1] >= inMon[i]) return false;
  }
  return true;
}


// deletes double variables, sorts
void normalizeMonom(monom_t &inMon) {
  // special cases 2
  if (inMon.s[0] > inMon.s[1]) { 
    uint32 tmp = inMon.s[0]; inMon.s[0] = inMon.s[1]; inMon.s[1] = tmp; // swap them
  } else 
    // same variable, x^2 = x -> insert one "1".
    if (inMon.s[0] == inMon.s[1]) { inMon.s[0] = 0; }  
  assert(isSortedMonom(inMon));
}


// returns the degree of a monomial
// in particular, ignores "1"-variables in a monomial
int monDeg(const monom_t &inMon) {
  assert(isSortedMonom(inMon));
  if (inMon[1] == 0) return 0;  // we cannot have degree -1 here (in this case, we do not have a monomial)
  if (inMon[0] == 0) return 1;
  return 2;
}


// comparison operators
bool operator < (const monom_t &lhs, const monom_t &rhs) {
  // only define the case maxPolyDeg == 2 
  assert((lhs[0] < lhs[1]) or ((lhs[0] == 0) and (lhs[1] == 0))); 
  assert((rhs[0] < rhs[1]) or ((rhs[0] == 0) and (rhs[1] == 0))); 

  // take care of the degree
  int leftDeg = monDeg(lhs); int rightDeg = monDeg(rhs);
  if (leftDeg != rightDeg) return leftDeg < rightDeg;

  // compare them
  if (lhs[1] < rhs[1]) return true;
  if (lhs[1] > rhs[1]) return false;   
  if (lhs[0] < rhs[0]) return true;
  if (lhs[0] >= rhs[0]) return false;   
  return false; // keep the compiler happy - we cannot reach this point actually...
}

bool operator > (const monom_t &lhs, const monom_t &rhs) {
  return rhs < lhs;  // re-use the operator from above
}

bool operator == (const monom_t &lhs, const monom_t &rhs) {
  // only define the case maxPolyDeg == 2
  return ((lhs.s[0] == rhs.s[0]) and (lhs.s[1] == rhs.s[1]));
}


// compare two monomials (sort-algorithm)
bool cmpMons (const monom_t &lhs, const monom_t &rhs) {
  return lhs > rhs;  // use operator from above
}


// again - but this time as a class (for priority_queue)
struct cmpMonsClass {
  bool operator() (const monom_t& lhs, const monom_t& rhs) const {
    // take care of the degree
    int leftDeg = monDeg(lhs); int rightDeg = monDeg(rhs);
    if (leftDeg != rightDeg) return leftDeg > rightDeg;

    if (lhs[1] > rhs[1]) return false;
    if (lhs[1] < rhs[1]) return true;
    if (lhs[0] > rhs[0]) return false;
    if (lhs[0] < rhs[0]) return true;
    return false; // they are equal here
  }
};


// hash operator for unordered_set & _map
struct hash_monom_t {
  size_t operator() (const monom_t& value) const {
    return value[0] ^ ((value[1] << 13) + (value[1] >> 13));
  }
};


bool isSortedMons(monVecQuad_t const &dest) {
  if (dest.size() == 0) return true;
  // check if all monomials are sorted
  for (size_t i = 0; i < dest.size(); i++) 
    if (!isSortedMonom(dest[i])) return false;

  if (dest.size() == 1) return true;
  // check if all monomials are in the right order
  // !!! use the comparison operator for monomaials here !!!
  for (size_t i = 1; i < dest.size(); i++) 
    if (!(dest[i-1] > dest[i])) return false;

  // O.K. - everything seems fine
  return true;
}


// copies the (normalized) version of src to dest. Overwrites "dest" in the process
void normalizeMons(monVecQuad_t &dest, monVecQuad_t &src) {
  assert(isSortedMons(src));
  dest.clear();  // delete dest
  if (src.size() == 0) return; // nothing to do
  // normalize all monomials
  for (size_t i = 0; i < src.size(); i++) normalizeMonom(src[i]);
  sort(src.begin(), src.end(), cmpMons);  // get src in shape
  size_t curPos = 0;
  // delete doubles in src
  while (curPos < src.size()-1) {
    if (src[curPos] == src[curPos+1]) {
      // two doubles - delete!
      curPos += 2;
    } else {
      // no doubles - we can copy an element to workingPoly
      dest.push_back(src[curPos]);
      curPos++;
    }
  }
  // make sure we do not forget the end
  if (curPos < src.size()) dest.push_back(src[curPos]);
}

// as above, but normalizes in place. Slightly less efficient, but avoids copying
void normalizeMons(monVecQuad_t &dest) {
  // normalize all monomials
  for (size_t i = 0; i < dest.size(); i++) normalizeMonom(dest[i]);
  if (dest.size() <= 1) return; // nothing to do
  // sort dest accordingly
  sort(dest.begin(), dest.end(), cmpMons);  
  size_t readPos = 0;
  size_t writePos = 0;
  // delete doubles in src
  while (readPos < dest.size()-1) {
    if (dest[readPos] == dest[readPos+1]) {
      // two doubles - delete!
      readPos += 2;
    } else {
      // no doubles - we can copy an element to workingPoly
      dest[writePos] = dest[readPos];
      readPos++; writePos++;
    }
  }
  // make sure we do not forget the end
  if (readPos < dest.size()) { dest[writePos] = dest[readPos]; writePos++; }
  
  // delete the last (readPos-writePos) elements
  if (readPos > writePos) dest.resize(writePos);
}

////////////////////////////////////////////////////////////////
///////////// varDict //////////////////////////////////////////
////////////////////////////////////////////////////////////////
// gives a mapping from variables (strings) to numbers (1..maxVarNum)
// !!! 1 is mapped to 0 !!!

const size_t maxVarNum = INT_MAX; // maximal number of variables

class varDict {
 public:
  // two maps to keep vars and numbers in sync
  map<string, uint32> var2num;
  map<uint32, string> num2var;
  // methods
  varDict() { clear();};
  void clear();
  uint32 giveNum(string var);
  string giveVar(uint32 num);
  size_t varsNum();
};

// constructor - assigns values to important variables
void varDict::clear() {
  var2num.clear();
  num2var.clear();

  // assign the "1"
  string workingString;
  workingString = "1";
  giveNum(workingString); // make sure this is a 0 for "1"
  assert(giveNum(workingString) == 0); // check
}

// get a number for each variable in question
// "1" -> 0
uint32 varDict::giveNum(string var) {
  // check if a variable already has a number
  // if not: assigns the first free number
  if (var2num.find(var) == var2num.end()) {
    assert(var2num.size() < maxVarNum);  // maximal number of variables in the system
    uint32 newNum = var2num.size(); 
    var2num[var] = newNum;
    num2var[newNum] = var;  // also think about the reverse direction
  }
  // O.K. - we have a number now. Return it!
  return var2num[var];
}

// get a variable for a given number - as long as it's defined
string varDict::giveVar(uint32 num) {
  assert((num >= 0) && (num < num2var.size()));
  return num2var[num];
}

// get the number of variables in the dictionary
size_t varDict::varsNum() {
  assert(var2num.size() != 0); 
  return var2num.size()-1;  // 1 is not strictly a variable...
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
// Define a variable
varDict fullVarDict;

//////////////////////////////////////////////////////////////////////////////////////////////////////
class linPoly {
 public:
  monVecType mons; 
  linPoly(); 
  linPoly(string _input);

  void clear() {mons.clear();}

  void setVal(uint32 *src); 

  bool isZero() const;
  bool isOne() const;
  bool isConstant() const;  

  bool isSorted() const;
  bool isEqual(linPoly &otherPoly) const;

  size_t len() const;
  uint32 lt() const;
  bool contains(uint32 tstMon) const;

  void add(linPoly &from);

  void printIt();
}; 

linPoly::linPoly() { clear(); }

// string constructor
linPoly::linPoly(string input) {
  clear();  // clear and nice object
  // parse the input, monomial for monomial
  bool addVar = false; 
  unsigned int pos = 0; string curVar(""); 
  // get all variables 
  while (pos < input.size()) {
    // extend current variable. Only allow a-z,A-Z,0-9 and '_'
    if (((input[pos] >= 'A') && (input[pos] <= 'Z')) || ((input[pos] >= 'a') && (input[pos] <= 'z')) ||
        ((input[pos] >= '0') && (input[pos] <= '9')) || (input[pos] == '_')) {curVar += input[pos];}
    if (input[pos] == '+') {addVar = true;}
    if (pos >= input.size()-1) addVar = true;
    // end of a variable
    if (addVar) {
      // do not add a zero
      if (curVar != "0") { 
        uint32 varNum = fullVarDict.giveNum(curVar);
        mons.push_back(varNum);
      }
      curVar = "";
      addVar = false; // need to wait for the next '+' sign
    }
    // house keeping
    pos++;
  }
  // normalize our result
  normalizeMons(mons);
}


void linPoly::setVal(uint32 *src) {
  // get the length
  uint32 inLen = src[0];
  mons.clear();
  for (uint32 i = 1; i <= inLen; i++) 
    mons.push_back(src[i]);  
  normalizeMons(mons);
}

bool linPoly::isZero() const {
  return this->len() == 0; 
}

bool linPoly::isOne() const {
  if (this->len() != 1) return false;
  return this->mons[0] == 0;
}

bool linPoly::isConstant() const {
  if (this->len() > 1) return false;
  if (this->len() == 0) return true;  // it's the zero polynomial
  return this->mons[0] == 0; // it's the one polynomial
}

bool linPoly::isSorted() const {
  if (this->len() <= 1) return true;
  for (size_t i= 0; i < this->len()-1; i++) 
    if (mons[i] < mons[i+1]) return false;
  return true;
}

bool linPoly::isEqual(linPoly &otherPoly) const {
  assert(isSorted()); 
  assert(otherPoly.isSorted());
  if (mons.size() != otherPoly.mons.size()) return false;
  for (size_t i = 0; i < mons.size(); i++)
    if (mons[i] != otherPoly.mons[i]) return false;
  return true;
}

size_t linPoly::len() const {
  return this->mons.size();
}

// returns the leading term of the polynomial
uint32 linPoly::lt() const {
  assert(len() != 0);
  return mons[0];
}

// checks if tstMon is containted in the current polynomial
bool linPoly::contains(uint tstMon) const {
  assert(isSorted());
  // use binary search
  int leftPos, rightPos, midPos;
  leftPos = 0;
  rightPos = len()-1;
  while (leftPos <= rightPos) {
    midPos = (leftPos+rightPos) / 2; 
    if (mons[midPos] == tstMon) return true;
    if (mons[midPos] > tstMon) 
      leftPos = midPos+1;
    else rightPos = midPos-1;
  }
  return false;
}

// adds two polynomials. Uses intermediate memory
void linPoly::add(linPoly &from) {
  assert(this->isSorted());
  assert(from.isSorted());
  monVecType destMons;
  size_t thisPos = 0, fromPos = 0;
  while ((thisPos < this->len()) and (fromPos < from.len())) {
    // skip as many monomials as possible on either side
    // copy monomials from "this" that are clearly before "from"
    while ((this->mons[thisPos] > from.mons[fromPos]) and (thisPos < this->len())) {
      destMons.push_back(this->mons[thisPos]);
      thisPos++;
    }
    if (thisPos >= this->len()) break;  // nothing left to do
    // copy monomials from "from" that are clearly before "this"
    while ((from.mons[fromPos] > this->mons[thisPos]) and (fromPos < from.len())) {
      destMons.push_back(from.mons[fromPos]);
      fromPos++;
    }
    if (fromPos >= from.len()) break; // nothing left to do
    // erase all double monomials ("1+1 = 0")
    while ((this->mons[thisPos] == from.mons[fromPos]) && (thisPos < this->len()) && (fromPos < from.len())) { 
      thisPos++; fromPos++; }
  }
  // copy the remaining monomials
  while (thisPos < this->len()) {
    destMons.push_back(this->mons[thisPos]); thisPos++; }
  while (fromPos < from.len()) {
    destMons.push_back(from.mons[fromPos]); fromPos++; }
  // copy all data back to "this"
  this->mons.swap(destMons);
}

void linPoly::printIt() {
  if (isZero()) cout << "0";
  for (size_t i = 0; i < this->len(); i++) {
    cout << fullVarDict.giveVar(this->mons[i]);
    //cout << " (" << this->mons[i] << ")";
    if (i != this->len()-1) cout << " + ";
  }
  cout << endl;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////
typedef unordered_map<uint32, linPoly> linPolyMat_t;
typedef unordered_set<uint32> monSet;
typedef unordered_map<uint32, monSet > polyByMon_t;

class linStore {
 public:
// variables
  linPolyMat_t linPolyMat;
  polyByMon_t polyByMon;
  bool isConsistent; // did we find an error in the overall system?

// constructor & near-constructor
  linStore() { clear(); }
  void clear();

// methods
  void insertPoly(linPoly & addPoly);
  size_t rank();
  size_t monNum();

  int giveValue(uint32 lt);
  int giveValue(string varName);

  bool checkConsistency();

  void printIt();
};

// deletes all internal variables
void linStore::clear() {
  linPolyMat.clear();
  polyByMon.clear();
  isConsistent = true;
}


// reduces and (if necessary) adds a polynomial to the store
void linStore::insertPoly(linPoly &addPoly) {
  // get some silly cases out of the way
  if (addPoly.isZero()) return;
  if (addPoly.isOne()) isConsistent = false;  // system is not consistent

  assert(checkConsistency()); // be prudent

  // reduce the polynomial by all monomials possible
  monVecType accuList;
  accuList.insert(accuList.end(),addPoly.mons.begin(), addPoly.mons.end());
  // do not check the lt; we do this below (after normalizing)
  for (size_t i = 1; i < addPoly.len(); i++) {
    uint32 curMon = addPoly.mons[i];
    // check if we can reduce by this monomial
    if (linPolyMat.find(curMon) != linPolyMat.end()) {
      // seems like -> copy!
      accuList.insert(accuList.end(), linPolyMat[curMon].mons.begin(), linPolyMat[curMon].mons.end());
    }
  }

  // normalize it (up to the leading term everything can have changed)
  linPoly workingPoly;
  normalizeMons(workingPoly.mons, accuList);

  // check if we are shorter than the current polynomial for the leading term
  // if yes: exchange them
  uint32 lt = workingPoly.lt(); 
  if (linPolyMat.find(lt) != linPolyMat.end()) {
    // seem this lt is already known. Reduce
    if (linPolyMat[lt].len() > workingPoly.len()) {
      // erase all terms from the corresponding sets (they will change)
      for (size_t i = 1; i < linPolyMat[lt].len(); i++) {
        uint32 curMon = linPolyMat[lt].mons[i];
        if (curMon == 0) break; // we are done
        polyByMon[curMon].erase(lt);  // mark that our polynomial contains this term
      }
      // we found a shorter equation => replace! 
      linPolyMat[lt].add(workingPoly);
      linPolyMat[lt].mons.swap(workingPoly.mons);  // swap the two vectors => swap the polynomials
      // add the terms we have new now
      for (size_t i = 1; i < linPolyMat[lt].len(); i++) {
        uint32 curMon = linPolyMat[lt].mons[i];
        if (curMon == 0) break; // we are done
        polyByMon[curMon].insert(lt);  // mark that our polynomial contains this term
      }
    } else workingPoly.add(linPolyMat[lt]);  // no - just reduce by this equation    
  }

  // did we generate a constant polynomial? 
  if (workingPoly.isZero()) return;  // we are done here
  if (workingPoly.isOne()) {
    isConsistent = false; // actually, not a real good outcome, complain!
    return;
  }

  // we have found a new basis element for the linear system
  // get the correct leading term
  lt = workingPoly.lt(); 
  linPolyMat[lt] = workingPoly;

  // register all monomials (except the first and "+1") in the corresponding set
  for (size_t i = 1; i < workingPoly.len(); i++) {
    uint32 curMon = workingPoly.mons[i];
    if (curMon == 0) break; // we are done
    polyByMon[curMon].insert(lt);  // mark that our polynomial contains this term
  }

  // reduce all other polynomials by this polynomial (if applicable)
  unordered_set< uint32 >::iterator pos;
  for (pos = polyByMon[lt].begin(); pos != polyByMon[lt].end(); ++pos) {
    uint32 otherLT = *pos;
    assert(otherLT != lt);
    // erase all terms from the corresponding sets
    for (size_t i = 1; i < linPolyMat[otherLT].len(); i++) {
      uint32 otherMon = linPolyMat[otherLT].mons[i];
      if (otherMon == 0) break; // we are done
      if (otherMon == lt) continue; // skip monomial that is skipped anyway
      polyByMon[otherMon].erase(otherLT);  // mark that our polynomial contains this term
    }
    // perform addition
    linPolyMat[otherLT].add(workingPoly);
    // insert all terms into the corresponding sets
    for (size_t i = 1; i < linPolyMat[otherLT].len(); i++) {
      uint32 otherMon = linPolyMat[otherLT].mons[i];
      if (otherMon == 0) break; // we are done
      polyByMon[otherMon].insert(otherLT);  // mark that our polynomial contains this term
    }
  }
  // empty set for leading term
  polyByMon.erase(lt);

  assert(checkConsistency());
}

// returns the rank of the store
size_t linStore::rank() {
  return linPolyMat.size();
}

// counts the total number of monomials in the store
// workload: linear in the number of polynomials, i.e. rank of the matrix
size_t linStore::monNum() {
  size_t res = 0;
  linPolyMat_t::const_iterator itr;
  for(itr = linPolyMat.begin(); itr != linPolyMat.end(); ++itr) {
    linPoly ourPoly = (*itr).second;
    res += ourPoly.len(); // get the number of monomials
  }
  return res;
}

// returns the value for a given leading term
// 0/1: zero/one constant
// 100: linear polynomial with more than one variable
//  -1: unknown
const int giveValue_unknown = -1;
const int giveValue_one = 1;
const int giveValue_zero = 0;
const int giveValue_linear = 100;
int linStore::giveValue(uint32 lt) {
  // check if we know this lt
  linPolyMat_t::const_iterator got = linPolyMat.find(lt);
  if (got == linPolyMat.end()) return giveValue_unknown;

  // make sure we do not have an invalid equation
  assert(linPolyMat[lt].lt() != 0);  
  assert(linPolyMat[lt].len() != 0);

  // find out in which case we are
  if ((linPolyMat[lt].len() == 1) && (linPolyMat[lt].lt() == lt)) return giveValue_zero;
  if ((linPolyMat[lt].len() == 2) && (linPolyMat[lt].mons[1] == 0)) return giveValue_one;
  return giveValue_linear;
}

// as above, but use a variable name
int linStore::giveValue(string varName) {
  return giveValue(fullVarDict.giveNum(varName));
}

// check if all internal data structures are 
// consistent with each other
bool linStore::checkConsistency() {

  // verify that all polynomials have the same lt as their key
  linPolyMat_t::const_iterator itr;
  for(itr = linPolyMat.begin(); itr != linPolyMat.end(); ++itr) {
    uint32 lt = (*itr).first;
    linPoly curPoly = (*itr).second;
    if (lt != curPoly.lt()) {
      cout << "Error with lt for " << lt << " / " << fullVarDict.giveVar(lt) << endl;
      cout << "  "; curPoly.printIt();
    }
  }
  // verify that all polynomials are fully reduced
  for(itr = linPolyMat.begin(); itr != linPolyMat.end(); ++itr) {
    linPoly curPoly = (*itr).second;
    for (uint32 i = 1; i < curPoly.len(); i++) {
      linPolyMat_t::const_iterator got = linPolyMat.find(curPoly.mons[i]);
      if (got != linPolyMat.end()) { cout << "Error with reduction "  << curPoly.mons[i] << " for "; curPoly.printIt(); }
    }
  }
  // verify that all monomials of each polynomial are registered in the corresponding set
  for(itr = linPolyMat.begin(); itr != linPolyMat.end(); ++itr) {
    linPoly curPoly = (*itr).second;
    uint32 lt = curPoly.lt();
    for (uint32 i = 1; i < curPoly.len(); i++) {
      uint32 curMon = curPoly.mons[i];
      if (curMon == 0) break; // we are done here
      unordered_set< uint32 >::iterator got = polyByMon[curMon].find(lt);
      if (got == polyByMon[curMon].end()) { 
        cout << "Error with registration "  << curPoly.mons[i] << " for "; curPoly.printIt(); 
      }
    }
  }
  // verify that all elements registered in sets are actually monomials in the corresponding polynomial
  polyByMon_t::const_iterator setItr;
  for (setItr = polyByMon.begin(); setItr != polyByMon.end(); ++setItr) {
    uint32 curMon = (*setItr).first;
    unordered_set< uint32 > curSet = (*setItr).second;
    unordered_set< uint32 >::const_iterator elemItr;
    for(elemItr = curSet.begin(); elemItr != curSet.end(); ++elemItr) {
      uint32 curLT = *elemItr;
      if (not(linPolyMat[curLT].contains(curMon))) { cout << "Found wrong element " << curMon << " for "; linPolyMat[curLT].printIt(); }
    }
  }
  return true;
}


void linStore::printIt() {
  cout << "rank: " << rank() << " / mons: " << monNum() << endl;
  linPolyMat_t::const_iterator itr;
  for(itr = linPolyMat.begin(); itr != linPolyMat.end(); ++itr) {
    linPoly ourPoly = (*itr).second;
    ourPoly.printIt(); 
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
linStore replaceLin;  // global variable for all linear equations


//////////////////////////////////////////////////////////////////////////////////////////////////////
// sparse polynomial over GF(2)
class polyGF2 {
 public:
  monVecQuad_t mons;

  // Constructors
  polyGF2() {};
  polyGF2(string inPoly);
  
  // make a new polynomial
  void clear();

  // some small, handy functions
  bool isZero() const { return mons.size() == 0; };
  size_t len() const { return mons.size(); }
  bool isSorted() const { return isSortedMons(mons); } 
  int deg() const; 
  monom_t lt() {assert(not(isZero())); return mons[0]; };

  // core functions
  int commonTerms (const polyGF2 &from) const;
  void add(monom_t &from);
  void add(polyGF2 const &otherPoly);
  bool eval(unsigned int maxOut, unsigned int maxIn); 
  bool sparseEval(unordered_set<monom_t, hash_monom_t> &monSet);

  // for debugging / output
  void printIt() const;
};

// string constructor
polyGF2::polyGF2(string _input) {
  clear(); // we want an empty polynomial
  if ((_input.length() == 0) || (_input == "0")) return; // nothing to do
  // parse the input, monomial for monomial
  // skip all space in the input (not efficient but works)
  string input(_input);
  unsigned int i = 0;
  while (i < input.size()) {
    if (input[i] == ' ') input.erase(i,1);
    else i++;
  } 
  // get everything ready for parsing
  bool addVar = false, addMon = false, isZero = false;
  unsigned int pos = 0; string curVar(""); 
  vector<uint32> curMon;
  // get all variables 
  while (pos < input.size()) {
    // extend current variable
    if (((input[pos] >= 'A') && (input[pos] <= 'Z')) || ((input[pos] >= 'a') && (input[pos] <= 'z')) ||
        ((input[pos] >= '0') && (input[pos] <= '9')) || (input[pos] == '_')) {curVar += input[pos];}
    if (input[pos] == '*') {addVar = true;}
    if (input[pos] == '+') {addVar = true; addMon = true;}
    if (pos >= input.size()-1) {addVar = true; addMon = true;}
    // end of a variable
    if (addVar) {
      // mark zero
      if (curVar == "0") isZero = true;  // there was a zero in the monomial
      else {
        uint32 varNum = fullVarDict.giveNum(curVar);
        curMon.push_back(varNum);
      }
      curVar = "";
      addVar = false; // need to wait for the next '*' sign
    }
    // end of monomial
    if (addMon) {
      // check if the degree is O.K.
      if (curMon.size() > maxPolyDeg) { cout << "Error - monomial too large at " << pos << " with input " << _input << endl; exit(1); }
      monom_t workingMon;
      for (size_t i = 0; i < maxPolyDeg; i++) {
        if (i < curMon.size()) workingMon[i] = curMon[i];
        else workingMon[i] = 0; // set all other variables to 0
      }
      if ((curMon.size() != 0) and not(isZero)) {     // only add non-zero monomials
        normalizeMonom(workingMon);
        mons.push_back(workingMon); 
      }
      curMon.clear(); // and now - space for a new monomial :-)
      addMon = false;  // wait for the next '+' sign
      isZero = false;
    }
    // house keeping
    pos++;
  }
  // normalize our result
  normalizeMons(mons);
}


bool operator == (const polyGF2 &lhs, const polyGF2 &rhs) {
  return lhs.mons == rhs.mons;
}


void polyGF2::clear() {
  mons.clear();
}


// -1: zero
//  0: 1 (constant)
//  1: deg 1
//  2: deg 2
int polyGF2::deg() const {
  assert(isSorted());
  if (isZero()) return -1;
  return monDeg(mons[0]);
}


// computes the number of common terms in this & from
int polyGF2::commonTerms(const polyGF2 &from) const {
  assert(this->isSorted());
  assert(from.isSorted());

  int resVal = 0;

  size_t thisPos = 0, fromPos = 0;
  while ((thisPos < this->len()) and (fromPos < from.len())) {
    // skip as many polynomials as possible on either side
    // skip monomials from "this" that are clearly before "from"
    while (this->mons[thisPos] < from.mons[fromPos]) {
      thisPos++;
      if (thisPos >= this->len()) return resVal;
    }
    if (thisPos >= this->len()) return resVal;  // nothing left to do
    // skip monomials from "from" that are clearly before "this"
    while (from.mons[fromPos] < this->mons[thisPos]) {
      fromPos++;
      if (fromPos >= from.len()) return resVal;
    }
    if (fromPos >= from.len()) return resVal; // nothing left to do
    // count all double monomials ("1+1 = 0")
    while ((this->mons[thisPos] == from.mons[fromPos]) && (thisPos < this->len()) && (fromPos < from.len())) { 
      thisPos++; fromPos++; resVal++; 
      if (thisPos >= this->len()) return resVal;
      if (fromPos >= from.len()) return resVal;
    }
  }
  return resVal;
}

// add the /from/ monomial to this (none-destructive to "from")
void polyGF2::add(monom_t &from) {
  assert(this->isSorted());
  monVecQuad_t destMons; 
  destMons.reserve(this->len() + 1); // avoid reallocation of memory
  // loop over all monomials in this vector
  size_t thisPos = 0;
  while (thisPos < this->len()) {
    // skip as many polynomials as possible

    // copy monomials from "this" that are clearly before "from"
    while (this->mons[thisPos] > from) {
      destMons.push_back(this->mons[thisPos]);
      thisPos++;
      if (thisPos >= this->len()) break;
    }
    if (thisPos >= this->len()) {
      destMons.push_back(from);
      break;
    }  // nothing left to do

    // copy monomials from "from" that are clearly before "this"
    if (from > this->mons[thisPos]) {
      destMons.push_back(from);
      break;
    }
    // skip if double monomial ("1+1 = 0")
    else { 
      thisPos++;
      break;
    }
  }

  // copy the remaining monomials
  while (thisPos < this->len()) {
    destMons.push_back(this->mons[thisPos]); thisPos++; }

  // copy all data back
  this->mons.swap(destMons);
}

// add the /from/ polynomial to this (none-destructive to "from")
void polyGF2::add(polyGF2 const &from) {
  assert(this->isSorted());
  assert(from.isSorted());
  monVecQuad_t destMons; 
  destMons.reserve(this->len() + from.len()); // avoid reallocation of memory
  // loop over all monomials in both vectors
  size_t thisPos = 0, fromPos = 0;
  while ((thisPos < this->len()) and (fromPos < from.len())) {
    // skip as many polynomials as possible on either side

    // copy monomials from "this" that are clearly before "from"
    while (this->mons[thisPos] > from.mons[fromPos]) {
      destMons.push_back(this->mons[thisPos]);
      thisPos++;
      if (thisPos >= this->len()) break;
    }
    if (thisPos >= this->len()) break;  // nothing left to do

    // copy monomials from "from" that are clearly before "this"
    while (from.mons[fromPos] > this->mons[thisPos]) {
      destMons.push_back(from.mons[fromPos]);
      fromPos++;
      if (fromPos >= from.len()) break;
    }
    if (fromPos >= from.len()) break; // nothing left to do

    // skip all double monomials ("1+1 = 0")
    while (this->mons[thisPos] == from.mons[fromPos]) { 
      thisPos++; fromPos++; 
      if ((thisPos >= this->len()) or (fromPos >= from.len())) break; // no data anymore
    }
  }

  // copy the remaining monomials
  while (thisPos < this->len()) {
    destMons.push_back(this->mons[thisPos]); thisPos++; }
  while (fromPos < from.len()) {
    destMons.push_back(from.mons[fromPos]); fromPos++; }

  // copy all data back
  this->mons.swap(destMons);
}

// multiplies two linear polynomials into a quadratic polynomial
// * first, second: source
// * target += first * second (add the result)
// the result is *not* normalized (needs to be done in another step)
// in addition: the leading term (=position 0) of first/second is ignored
// note: result needs to be normalized outside
// toDo: make this more generic for degree 3 & 4
void mulAddLin(monVecQuad_t &target, linPoly &first, linPoly &second) {
  assert(maxPolyDeg == 2); // code below only works for deg=2
  for (size_t i = 1; i < first.len(); i++) {
    for (size_t j = 1; j < second.len(); j++) {
      monom_t curMon;
      curMon[0] = first.mons[i];
      curMon[1] = second.mons[j];
      normalizeMonom(curMon); // ensure right sorting     
      target.push_back(curMon);
    }
  }
}

// transfers a linear polynomial into a quadratic one.
// first term ("lt") is ignored
void lin2higher(monVecQuad_t &target, linPoly &from, uint32 otherVar = 0) {
  assert(maxPolyDeg == 2); // code only works for degree 2
  monom_t curMon;
  for (size_t i = 1; i < from.len(); i++) {
    curMon[0] = otherVar; // mostly 1-factor
    curMon[1] = from.mons[i];
    normalizeMonom(curMon); // ensure right sorting
    target.push_back(curMon);
  }
}

// copies a quadratic polynomial into a linear one
// assumes that the target polynomial
// * is normalized
// * has degree 1
// includes the leading term
// overwrites "target"
void higher2lin(linPoly &target, monVecQuad_t &from) {
  assert(isSortedMons(from));
  target.clear();  // free some space
  if (from.size() == 0) return;  // nothing to do
  assert(monDeg(from[0]) <= 1); // just to be sure

  // copy all variables into the other data structure
  for (size_t i = 0; i < from.size(); i++) {
    assert(from[i][0] == 0);  // check if it is really degree 1
    target.mons.push_back(from[i][1]);
  } 
}

// evaluates the given polynomial using linear equations from 
// current polynomial is replaced / destroyed
// return true: we have changed something
// * outMax: maximal number of newly added monomials
// * inMax: maximal size for the input polynomial (quadratic case)
bool polyGF2::eval(unsigned int outMax, unsigned int inMax) {
  monVecQuad_t accu; 
  bool retVal = false;
  size_t pos = 0;
  unsigned int addMons = 0;

  // start with degree 2 monomials
  while ((pos < len()) and (monDeg(mons[pos]) == 2)) {
    // go through all 4 possible cases for known/unknown
    bool knownZero = (replaceLin.giveValue(mons[pos][0]) != giveValue_unknown);
    bool knownOne = (replaceLin.giveValue(mons[pos][1]) != giveValue_unknown);
    retVal |= knownZero | knownOne;  // did we change anything?
    // both are known -> multiply (except they are too long)
    if (knownZero and knownOne) {
      // first, check we are allowed to evaluate
      if (replaceLin.linPolyMat[mons[pos][0]].len() >= inMax) return false;
      if (replaceLin.linPolyMat[mons[pos][1]].len() >= inMax) return false;
      addMons += replaceLin.linPolyMat[mons[pos][0]].len() 
               * replaceLin.linPolyMat[mons[pos][1]].len();
      if (addMons >= outMax) return false;
      // yes - we are :-)
      mulAddLin(accu, replaceLin.linPolyMat[mons[pos][0]], replaceLin.linPolyMat[mons[pos][1]] );
    }
    // none is known -> just copy
    else if (not(knownZero) and not(knownOne)) 
      accu.push_back(mons[pos]);
    else {
    // only one is known -> select accordingly
      uint32 knownVar, unknownVar;
      if (knownOne) { knownVar = mons[pos][1]; unknownVar = mons[pos][0]; }
      else {knownVar = mons[pos][0]; unknownVar = mons[pos][1]; }
      lin2higher(accu, replaceLin.linPolyMat[knownVar], unknownVar );
    }
    pos++; // done - we can deal with the next monomial now
  }
  // degree 1 monomials
  while ((pos < len()) and (monDeg(mons[pos]) == 1)) {
    assert(isSortedMonom(mons[pos]));
    assert(mons[pos][0] == 0);
    uint32 var = mons[pos][1];
    if (replaceLin.giveValue(var) == giveValue_unknown)
      accu.push_back(mons[pos]);  // no value known -> copy the current monomial!
    else {
        retVal = true;
	// we have an element in store - add to result
        lin2higher(accu, replaceLin.linPolyMat[var]);  
    }
    pos++; 
  }  
  // constant term (if any)
  if (pos < len())  {
    assert(monDeg(mons[pos]) == 0);  // we must have a *one* here
    accu.push_back(mons[pos]);  
    pos++;
  }
  assert(pos == len());
  if (accu.size() != 0) {
    for (size_t i = 0; i < accu.size(); i++) normalizeMonom(accu[i]);
    normalizeMons(mons, accu); 
  } else mons.clear();

  return retVal;
}

// evaluates the given polynomial using linear equations from 
// current polynomial is replaced / destroyed
// return true: we have changed something
bool polyGF2::sparseEval(unordered_set<monom_t, hash_monom_t> &monSet) {
  monVecQuad_t accu, fromList; 
  bool retVal = false;
  unordered_set<uint32> varSet; varSet.rehash(2*len());
  // copy the current polynomial
  for (auto monsIter = mons.begin(); monsIter != mons.end(); ++monsIter) {
    fromList.push_back(*monsIter);
    if ((*monsIter)[0] != 0) varSet.insert((*monsIter)[0]);
    if ((*monsIter)[1] != 0) varSet.insert((*monsIter)[1]);
  }
  vector<uint32> varList; 
  for (auto itr = varSet.begin(); itr != varSet.end(); ++itr) {
      // check if we know this variable
      if (replaceLin.giveValue(*itr) == giveValue_unknown) continue; 
      if (replaceLin.linPolyMat[*itr].len() > 2) continue;
    varList.push_back(*itr);
  }

  for (auto varIter = varList.begin(); varIter != varList.end(); ++varIter) {
    uint32 curVar = *varIter;
    accu.clear();
    size_t pos = 0;
    // start with degree 2 monomials
    while ((pos < len()) and (monDeg(mons[pos]) == 2)) {
      uint32 knownVar = 0; uint32 unknownVar = 0;
      assert(fromList[pos][0] != fromList[pos][1]);
      if (fromList[pos][0] == curVar) { 
        knownVar = curVar; unknownVar = fromList[pos][1]; 
      }
      if (fromList[pos][1] == curVar) { 
        knownVar = curVar; unknownVar = fromList[pos][0]; 
      }
      // only one is known -> select accordingly
      if (knownVar != 0) lin2higher(accu, replaceLin.linPolyMat[knownVar], unknownVar );
      else accu.push_back(mons[pos]);
      pos++; // done - we can deal with the next monomial now
    }
    // degree 1 monomials
    while ((pos < len()) and (monDeg(fromList[pos]) == 1)) {
      assert(isSortedMonom(fromList[pos]));
      assert(fromList[pos][0] == 0);
      uint32 var = fromList[pos][1];
      if (replaceLin.giveValue(var) == giveValue_unknown) {
        accu.push_back(fromList[pos]);  // no value known -> copy the current monomial!
      } else {
	// we have an element in store - add to result
        lin2higher(accu, replaceLin.linPolyMat[var]);  
      }
      pos++; 
    }  
    // constant term (if any)
    if (pos < len())  {
      assert(monDeg(fromList[pos]) == 0);  // we must have a *one* here
      accu.push_back(fromList[pos]);  
      pos++;
    }
    assert(pos == len());
    if (accu.size() != 0) {
      for (size_t i = 0; i < accu.size(); i++) normalizeMonom(accu[i]);
      normalizeMons(accu); 
      // is it a good polynomial?
      pos = 0; bool goodPoly = true;
      while ((pos < accu.size()) and goodPoly) {
        if (monDeg(accu[pos]) < 2) break;
	if (monSet.count(accu[pos]) == 0) goodPoly = false;
	pos++;
      }
      if (goodPoly) {
         retVal = true;
         mons.swap(accu);  // yes -> copy!
	 fromList.clear();
         for (auto monsIter = mons.begin(); monsIter != mons.end(); ++monsIter) 
           fromList.push_back(*monsIter);	 
      } else fromList.swap(accu);	 
    } else mons.clear();
  } // end for var
//  fromList.clear();
//  accu.clear();
  return retVal;
}


void polyGF2::printIt() const {
  if (isZero()) cout << "0";
  for (size_t i = 0; i < this->len(); i++) {
    if (i != 0) cout << " + ";
    // only output the part of the monomial that we need
    if (mons[i][0] != 0) {
      cout << fullVarDict.giveVar(this->mons[i][0]) << "*";
      cout << fullVarDict.giveVar(this->mons[i][1]); 
    } else cout << fullVarDict.giveVar(this->mons[i][1]);
    // also output the numbers
    //cout << " (" << mons[i][0] << ":" << mons[i][1] << ") ";    
  }
  cout << endl;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

typedef vector<polyGF2> quadPolyMat_t;
typedef unordered_set<uint32> idSet_t;
typedef unordered_map<monom_t, idSet_t, hash_monom_t > quadPolyByMon_t;
typedef unordered_map<uint32, idSet_t > quadPolyByVar_t;
typedef stack<int> freePos_t;


class quadStore_t {
 public:
// variables of the class
  quadPolyMat_t quadPolyList;
  quadPolyByMon_t quadPolyByMon;
  quadPolyByVar_t quadPolyByVar;
  freePos_t emptyRows;

// constructors and near-constructors
  quadStore_t() { clear(); };
  void clear();

// methods
  // internal methods --- only call when you know what you are doing
  int joinPoly(polyGF2 & inPoly, int hashPos);     // copies the full polynomial
  void removePoly(uint32 idx, bool fullRemove);        // works on its index

  bool commonVars(size_t startVal, size_t endVal, set<uint32> &resSet,
       monVecQuad_t workingList, unordered_set<monom_t, hash_monom_t> &monSet);
  bool halfCommonVars(size_t startVal, size_t endVal, set<uint32> &resSet,
       monVecQuad_t &workingList);
  

  // main interface for new polynomials
  // does a full reduction by other polynomials, too
  void insertPoly(polyGF2 & addPoly);
  void softReduce(size_t maxMonNum, size_t maxPolyInc);
  void echelonizeSparse(size_t maxMonNum);
  void echelonizeDense();
  size_t fullEval(size_t maxMonNum, size_t outMax, size_t inMax);
  size_t fullSparseEval();
  size_t sl();
  size_t hsl(float density, float densityMax);

  // statistics for the outside world
  size_t polyCnt() { assert(emptyRows.size() <= quadPolyList.size()); return quadPolyList.size() - emptyRows.size(); };
  size_t monNum();
  void replaceStat(uint32 &cntVal, uint32 &minVal, uint32 &maxVal, uint32 &medVal, float &avgVal);

  // make sure all data structures are in sync
  bool checkConsistency();

  void printIt(bool withNum);
};

// near-constructor: Brings everything back to normal, i.e. deletes all data structures
void quadStore_t::clear() {
  quadPolyList.clear();
  quadPolyByMon.clear();
  quadPolyByVar.clear();  
  // unfortunately, a stack does not have a "clear" method -> we do this ourselves
  for (size_t i = 0; i < emptyRows.size(); i++) 
    emptyRows.pop(); // just removes, do not return element
}


// adds a polynomial to the overall data structure
// does not reduce this polynomial
// returns: position of the polynomial in quadPolyList
int quadStore_t::joinPoly(polyGF2 & addPoly, int hasPos = -1) {
  // make sure the polynomial is in shape
  assert(addPoly.isSorted());

  // check initial conditions
  if (addPoly.isZero()) return -1;  // nothing to do
  // linear polynomial -> add to replaceLin!
  if (addPoly.deg() <= 1) {
    linPoly tmpLinPoly;
    higher2lin(tmpLinPoly, addPoly.mons);
    replaceLin.insertPoly(tmpLinPoly);
    return -1;
  }
  
  // get a new position for this polynomial (if needed)
  int myPos;
  if (hasPos >= 0) {
    myPos = hasPos;
    assert(addPoly == quadPolyList[myPos]);
  } else {
    // do we still have space in the vector? 
    if (emptyRows.size() != 0) {
      // there is an empty row -> use it!
      myPos = emptyRows.top();
      emptyRows.pop();
      quadPolyList[myPos] = addPoly;
    } else {
      // add this element to the (new) last position in the vector
      myPos = (int)quadPolyList.size();
      quadPolyList.push_back(addPoly);  // insert at this position (is a new polynomial)
      assert(quadPolyList[myPos] == addPoly);
    } 
  } 
  // now as we have a position in the list, 
  // get all sets up to date
  // register in monomial set (only deg2-mons)  
  for (size_t i = 0; i < addPoly.len(); i++) {
    monom_t curMon = addPoly.mons[i];
    if (monDeg(curMon) <= 1) break; // we are done
    quadPolyByMon[curMon].insert(myPos);  // mark that our polynomial contains this term
  }

  // register in var set
  unordered_set<uint32> varSet;
  varSet.rehash(addPoly.len()*2); // make enough space 
  for (size_t i = 0; i < addPoly.len(); i++) {
    monom_t curMon = addPoly.mons[i];
    // only add variables, not 1-factors
    if (curMon[0] != 0) varSet.insert(curMon[0]);
    if (curMon[1] != 0) varSet.insert(curMon[1]);
  }
  // add "myPos" to all variables
  for ( auto it = varSet.begin(); it != varSet.end(); ++it ) 
    quadPolyByVar[*it].insert(myPos);

  // we are done -> tell our position to the caller
  return myPos;
}


// frees the given index
// removes the index from polyByMon and polyByVar
//  * fullRemove: false: leave the polynomial in the list
void quadStore_t::removePoly(uint32 idx, bool fullRemove = true) {
  polyGF2 workingPoly = quadPolyList[idx];

  // is this a "real" polynomial?
  assert(not(workingPoly.isZero())); 
  assert(workingPoly.deg() == 2); 

  // unregister from "varSet"
  // first: collect all variables in a set
  unordered_set<uint32> varSet;
  varSet.rehash(workingPoly.len()*2); // make space 
  for (size_t i = 0; i < workingPoly.len(); i++) {
    monom_t curMon = workingPoly.mons[i];
    // only add variables, not 1-factors
    if (curMon[0] != 0) varSet.insert(curMon[0]);
    if (curMon[1] != 0) varSet.insert(curMon[1]);
  }
  // remove "myPos" for all variables
  for ( auto it = varSet.begin(); it != varSet.end(); ++it ) {
    uint32 varIdx = *it;
    quadPolyByVar[varIdx].erase(idx);
    if (quadPolyByVar[varIdx].size() == 0) quadPolyByVar.erase(varIdx);
  }

  // unregister in monomial set (only deg2-mons)
  for (size_t i = 0; i < workingPoly.len(); i++) {
    monom_t curMon = workingPoly.mons[i];
    if (monDeg(curMon) <= 1) break; // we are done
    quadPolyByMon[curMon].erase(idx);  // remove our polynomial from the set
    if (quadPolyByMon[curMon].size() == 0) quadPolyByMon.erase(curMon);
  }

  // clear the row, mark as empty
  if (fullRemove) { quadPolyList[idx].clear(); emptyRows.push(idx); }
}


// adds a polynomial to the overall data structure
// performs no reduction (see method echelonze)
void quadStore_t::insertPoly(polyGF2 &addPoly) {
  assert(checkConsistency());  // make sure everything is fine

  // check initial conditions
  if (addPoly.isZero()) return;  // nothing to do
  // linear polynomial -> add to replaceLin!
  if (addPoly.deg() <= 1) {
    linPoly tmpLinPoly;
    higher2lin(tmpLinPoly, addPoly.mons);
    replaceLin.insertPoly(tmpLinPoly);
    return;  // nothing to do for us
  }

  // it's degree 2, to put it into the matrix
  joinPoly(addPoly);  // that's all
}


// Performs echelonization on the polynomials in this
// stops if the number of monomials is higher than "maxMonNum"
void quadStore_t::echelonizeSparse(size_t maxMonNum = (1 << 30)) {
  assert(checkConsistency());  // make sure everything is fine

  size_t currentMonNum = monNum(); 
  // check if there is anything to do
  if (currentMonNum >= maxMonNum) return;

  // step 1: get the leading terms of all polynomials (only deg 2)
  // assumption: we have less polynomials than deg-2-monomials
  priority_queue<monom_t, deque<monom_t>, cmpMonsClass> largestLT;

  for (auto iter = quadPolyList.begin(); iter != quadPolyList.end(); ++iter) {
    if (iter->deg() <= 1) continue;
    monom_t lt = iter->lt();
    assert(monDeg(lt) >= 2); // better save than sorry...
    // only include lt that can be used for reduction
    if (quadPolyByMon[lt].size() <= 1) continue; 
    largestLT.push(lt);
  }

  // loop as long as we have some lt left
  while (!largestLT.empty()) {
    monom_t curLT;
    // get the largest element from the top --- we are non-empty here
    curLT = largestLT.top();
    assert(monDeg(curLT) == 2);
    // remove doubles in the priority queue
    while (!largestLT.empty()) {
      if (largestLT.top() == curLT) largestLT.pop();
      else break;
    }

    // take all polynomials that contain this monomial. Extract the ones that have it as a lt. Search for shortest element
    uint32 minLen = INT_MAX; uint32 minPolyNum = INT_MAX;
    for (auto iter = quadPolyByMon[curLT].begin(); iter != quadPolyByMon[curLT].end(); ++iter) {
      polyGF2 workingPoly = quadPolyList[*iter];
      assert(!workingPoly.isZero());
      if ((workingPoly.lt() == curLT) and (workingPoly.len() < minLen)) {
        minLen = workingPoly.len();
        minPolyNum = *iter;
      }
    }

    // hopefully, we have a minimum. If not - continue
    if (minLen == INT_MAX) continue;

    // copy all numbers for this monomial (will cause problems with the iterator otherwise)
    vector <uint32> workingList;
    workingList.reserve(quadPolyByMon[curLT].size());
    for (auto iter = quadPolyByMon[curLT].begin(); iter != quadPolyByMon[curLT].end(); ++iter) 
      if (*iter != minPolyNum) workingList.push_back(*iter);

    // reduce all other polynomials by this polynomial
    polyGF2 echelonElem = quadPolyList[minPolyNum];

    // go over the list of copied items, reduce all elements (except echelon element)
    for (auto iter = workingList.begin(); iter != workingList.end(); ++iter) {
      uint32 workingNum = (*iter);  
      assert(workingNum != minPolyNum);

      // delete external data structure
      currentMonNum -= quadPolyList[workingNum].len();
      removePoly(workingNum, false);

      // reduce by the echelon element
      quadPolyList[workingNum].add(echelonElem);

      // get some special cases out of the way
      if (quadPolyList[workingNum].deg() <= 1) {
        joinPoly(quadPolyList[workingNum], workingNum);  // use a side-effect here: will copy the polynomial to linStore (except: it's zero)
        quadPolyList[workingNum].clear();  // poly is gone -> make space
        emptyRows.push(workingNum);  // also register it as empty
        continue; // next poly
      }

      // record the new lt (if smaller than current one and deg >= 2)
      assert(monDeg(quadPolyList[workingNum].lt()) >= 2);
      if (curLT > quadPolyList[workingNum].lt()) largestLT.push(quadPolyList[workingNum].lt());

      // reintegrate the polynomial into the matrix
      joinPoly(quadPolyList[workingNum], workingNum);
      currentMonNum += quadPolyList[workingNum].len();
      // make sure we don't doo too much
      if (currentMonNum >= maxMonNum) return;
    }
  }
  assert(checkConsistency());  // check if everything is still fine


  /*
  // below is for debug only
  cout << "I am here" << endl; 

  // ensure that all reductions were performed (set look-up)
  for (auto iter = quadPolyList.begin(); iter != quadPolyList.end(); ++iter) {
    if (iter->deg() <= 1) continue; // get rid of zero 1 polynomials
    monom_t lt = iter->lt();
    assert(quadPolyByMon[lt].size() == 1); 
    if (quadPolyByMon[lt].size() != 1) {
      cout << "Error with set ";
      for (auto iter = quadPolyByMon[lt].begin(); iter != quadPolyByMon[lt].end(); ++iter)
        cout << *iter << ", ";
      cout << endl;
    }
  }

  // again - but this time monomial by monomial
  for (auto ltPolyIter = quadPolyList.begin(); ltPolyIter != quadPolyList.end(); ltPolyIter++) {
    if (ltPolyIter->deg() <= 1) continue; // get rid of zero 1 polynomials
    monom_t testLT = ltPolyIter->lt();
    for (auto polyIter = quadPolyList.begin(); polyIter != quadPolyList.end(); polyIter++) {
      if (polyIter->deg() <= 1) continue; // get rid of zero 1 polynomials
      for (auto monIter = polyIter->mons.begin()++; monIter != polyIter->mons.end(); monIter++) {
        if (((*monIter)[0] == 90) and ((*monIter)[1] == 107)) cout << "got poly: "; polyIter->printIt();
        if (testLT == *monIter) {
          cout << "found problem with " << (*monIter)[0] << ":" << (*monIter)[1] << " and" << endl;
          cout << "  first poly: "; polyIter->printIt();
          cout << "     lt poly: "; ltPolyIter->printIt();
        }
      }
    }
  }
  */
}


// Performs echelonization on the polynomials in this
// stops if the number of monomials is higher than "maxMonNum"
void quadStore_t::softReduce(size_t maxMonNum = (1 << 30), size_t maxPolyExpand = 0) {
  assert(checkConsistency());  // make sure everything is fine

  // check if there is anything to do
  size_t currentMonNum = monNum();
  if (currentMonNum >= maxMonNum) return;

  // loop over all polynomials
  for (size_t curPos = 0; curPos < quadPolyList.size(); curPos++) {
    if (quadPolyList[curPos].len() == 0) continue; // nothing to do here
    const polyGF2 reducerPoly = quadPolyList[curPos]; 
    // search for all polynomials that have at least one quadratic monomial in common
    set<uint32> otherPolys;
    for (auto iter = reducerPoly.mons.begin(); iter != reducerPoly.mons.end(); ++iter) {
      const monom_t targetMon = *iter;
      if (monDeg(targetMon) <= 1) break; // we are done here
//      int insCnt = 0;
      for (auto iterPoly = quadPolyByMon[targetMon].begin(); iterPoly != quadPolyByMon[targetMon].end(); ++iterPoly) {
        otherPolys.insert(*iterPoly); // insert all polynomials 
//        insCnt++;
//        if (insCnt > 10) break; // don't register too much
      }
    }

    // loop over all polynomials and try to reduce the number of quadratic polynomials
    for (auto polyNumIter = otherPolys.begin(); polyNumIter != otherPolys.end(); ++polyNumIter) {
      const size_t polyNum = *polyNumIter;
      if (polyNum == curPos) continue; // don't reduce yourself
      uint32 commonTerms = reducerPoly.commonTerms(quadPolyList[polyNum]);
      if (reducerPoly.len() >= (2*commonTerms+maxPolyExpand)) continue;  // result will be bigger -> return
      // some book-keeping
      currentMonNum -= quadPolyList[polyNum].len();
      removePoly(polyNum, false);
      // reduce the polynomial
      quadPolyList[polyNum].add(reducerPoly);
      // get some special cases out of the way
      if (quadPolyList[polyNum].deg() <= 1) {
        joinPoly(quadPolyList[polyNum], polyNum);  // use a side-effect here: will copy the polynomial to linStore (except: it's zero)
        quadPolyList[polyNum].clear();  // poly is gone -> make space
        emptyRows.push(polyNum);  // also register it as empty
        continue; // next poly
      }
      assert(monDeg(quadPolyList[polyNum].lt()) >= 2);
      // reintegrate the polynomial into the matrix
      joinPoly(quadPolyList[polyNum], polyNum);
      currentMonNum += quadPolyList[polyNum].len();
      // make sure we don't doo too much
      if (currentMonNum >= maxMonNum) return;
    }
  }
  assert(checkConsistency());  // check if everything is still fine
}


// uses M4RI echelonization. Mainly converts from one 
void quadStore_t::echelonizeDense() {
  assert(checkConsistency());  // check if everything is still fine
  // create a list of monomials
  vector<monom_t> semantic;
  semantic.reserve( quadPolyByMon.size() + quadPolyByVar.size() + 1); // make sure we have enough space
  // start with the quadratic parts of all polynomials
  for (auto curItem = quadPolyByMon.begin(); curItem != quadPolyByMon.end(); ++curItem) 
    semantic.push_back(curItem->first);  // add the current monomial
  // "+1"
  monom_t tmpMon; tmpMon[0] = 0; tmpMon[1] = 0;
  semantic.push_back(tmpMon);
  // continue with the linear parts of all polynomials (slight overkill - but won't hurt)
  for (auto curItem = quadPolyByVar.begin(); curItem != quadPolyByVar.end(); ++curItem) {
    tmpMon[1] = curItem->first;
    assert(isSortedMonom(tmpMon));
    semantic.push_back(tmpMon);  // add the current monomial
  }
  // sort according to our ordering
  sort(semantic.begin(), semantic.end(), cmpMons);  
  assert(isSortedMons(semantic));
  // create a reverse look-up table
  unordered_map<monom_t, uint32, hash_monom_t> reverseSemantic;
  reverseSemantic.rehash(semantic.size());
  for (uint32 curPos = 0; curPos < semantic.size(); curPos++) 
    reverseSemantic[semantic[curPos]] = curPos;

  // create a M4RI matrix of size (#poly, #semantic)
  mzd_t *denseMatrix = mzd_init(quadPolyList.size()-emptyRows.size(), semantic.size());
 
  // check: is this matrix empty?
  for (rci_t row = 0; row < denseMatrix->nrows; row++) 
    for (rci_t col = 0; col < denseMatrix->ncols; col++) 
      assert(mzd_read_bit(denseMatrix, row,col) == 0);

  // copy all non-zero entries into the M4RI-object
  uint32 writeRow = 0;
  for (uint32 readRow = 0; readRow < quadPolyList.size(); readRow++) {
    if (quadPolyList[readRow].len() == 0) continue;  // nothing to read
    for (auto curMon = quadPolyList[readRow].mons.begin(); curMon != quadPolyList[readRow].mons.end(); ++curMon) {
      assert(reverseSemantic.find(*curMon) != reverseSemantic.end());
      mzd_write_bit(denseMatrix, writeRow, reverseSemantic[*curMon], 1);
    }
    removePoly(readRow, true); // fully delete polynomial
    writeRow++;  // next row
  }

  // are all polynomials deleted?
  assert(emptyRows.size() == quadPolyList.size());

  // echelonize (1: full echelon form / 0: auto select size)
  rci_t rank = mzd_echelonize_m4ri(denseMatrix, 1, 0);

  // copy all data back into quadPolyList
  for (rci_t row = 0; row < rank; row++) {
    polyGF2 tmpPoly; tmpPoly.clear();
    for (rci_t col = 0; col < denseMatrix->ncols; col++) 
      if (mzd_read_bit(denseMatrix, row,col)) tmpPoly.mons.push_back(semantic[col]);    
    // add polynomial to list
    if (tmpPoly.len() != 0) {
      assert(tmpPoly.isSorted());
      joinPoly(tmpPoly);
    }
  }
  // free space
  mzd_free(denseMatrix);

  assert(checkConsistency());  // check if everything is still fine
}


// evaluate all polynomials in store
// returns the number of polynomials changed
// * maxMonNum: maximal size of the overall system
// * outMax: maximal number of monomials that can be added
// * inMax: maximal number of monomials that can be used as in input for multiplication
size_t quadStore_t::fullEval(size_t maxMonNum = (1 << 31), size_t outMax = (1 << 31) , size_t inMax = (1 << 31)) {
  size_t retNum = 0;

  assert(checkConsistency());  // make sure everything is fine

  // check if there is anything to do
  size_t currentMonNum = monNum();
  if (currentMonNum >= maxMonNum) return 0;

  // first: compute the intersection between variables in use and known variables
  // second: compute the polynomials we need to change
  idSet_t polySet; polySet.rehash(replaceLin.linPolyMat.size());
  if (replaceLin.rank() < quadPolyByVar.size()) {
    // assume: the number of linear equations is smaller
    for (auto itr = replaceLin.linPolyMat.begin(); itr != replaceLin.linPolyMat.end(); ++itr) {
      auto foundPos = quadPolyByVar.find(itr->first);
      if (foundPos == quadPolyByVar.end()) continue; // nothing to do
      // add all polynomials to polySet
      for (auto setItr = foundPos->second.begin(); setItr != foundPos->second.end(); ++setItr) 
        polySet.insert(*setItr);
    }
  } else {
    // assume: the number of variables in use is smaller
    for (auto itr = quadPolyByVar.begin(); itr != quadPolyByVar.end(); ++itr) {
      if (replaceLin.giveValue(itr->first) == giveValue_unknown) continue; // check if we know this variable
      // seems like => add polynomials to set
      for (auto setItr = itr->second.begin(); setItr != itr->second.end(); ++setItr) 
        polySet.insert(*setItr);
    }
  }

  // evaluate all polynomials in store
  for (auto setItr = polySet.begin(); setItr != polySet.end(); ++setItr) {
    uint32 const i = *setItr; // make the code below a bit more readable
    assert(not(quadPolyList[i].isZero())); // to be on the safe side
    // first: unregister the polynomial
    currentMonNum -= quadPolyList[i].len();
    removePoly(i, false);
    // second: fully evaluate it
    if (quadPolyList[i].eval(outMax, inMax)) retNum++;
    // get two special cases out of the way: 0 or lin
    if (quadPolyList[i].deg() <= 1) {
        joinPoly(quadPolyList[i], i);  // use a side-effect here: will copy the polynomial to linStore (except: it's zero)
        quadPolyList[i].clear();  // poly is gone -> make space
        emptyRows.push(i);  // also register it as empty
        continue; // next poly - nothing to do here
    }
    // third: re-add the "new" polynomial
    joinPoly(quadPolyList[i], i);
    // check if we are done here
    currentMonNum += quadPolyList[i].len();
    if (currentMonNum >= maxMonNum) return retNum;
  }

  assert(checkConsistency());  // make sure everything is still fine

  return retNum;
}

// evaluate all polynomials in store
// returns the number of polynomials changed
// * maxMonNum: maximal size of the overall system
// * outMax: maximal number of monomials that can be added
// * inMax: maximal number of monomials that can be used as in input for multiplication
size_t quadStore_t::fullSparseEval() {
  size_t retNum = 0;

  assert(checkConsistency());  // make sure everything is fine

  // generate a set of all allowed (quadratic) monomials
  unordered_set< monom_t, hash_monom_t > monSet;
  for (auto itr = quadPolyByMon.begin(); itr != quadPolyByMon.end(); ++itr)
    monSet.insert((*itr).first);
  for (auto itr = replaceLin.polyByMon.begin(); itr != replaceLin.polyByMon.end(); ++itr){
    monom_t curMon;
    curMon[0] = 0;
    curMon[1] = (*itr).first;
    monSet.insert(curMon);
  }

  // first: compute the intersection between variables in use and known variables
  // second: compute the polynomials we need to change
  idSet_t polySet; polySet.rehash(replaceLin.linPolyMat.size());
  if (replaceLin.rank() < quadPolyByVar.size()) {
    // assume: the number of linear equations is smaller
    for (auto itr = replaceLin.linPolyMat.begin(); itr != replaceLin.linPolyMat.end(); ++itr) {
      auto foundPos = quadPolyByVar.find(itr->first);
      if (foundPos == quadPolyByVar.end()) continue; // nothing to do
      // add all polynomials to polySet
      for (auto setItr = foundPos->second.begin(); setItr != foundPos->second.end(); ++setItr) 
        polySet.insert(*setItr);
    }
  } else {
    // assume: the number of variables in use is smaller
    for (auto itr = quadPolyByVar.begin(); itr != quadPolyByVar.end(); ++itr) {
      if (replaceLin.giveValue(itr->first) == giveValue_unknown) continue; // check if we know this variable
      // seems like => add polynomials to set
      for (auto setItr = itr->second.begin(); setItr != itr->second.end(); ++setItr) 
        polySet.insert(*setItr);
    }
  }
  
  // evaluate all polynomials in store
  for (auto setItr = polySet.begin(); setItr != polySet.end(); ++setItr) {
    uint32 const i = *setItr; // make the code below a bit more readable
    assert(not(quadPolyList[i].isZero())); // to be on the safe side
    // first: unregister the polynomial
    removePoly(i, false);
    // second: fully evaluate it
    if (quadPolyList[i].sparseEval(monSet)) retNum++;
    // get two special cases out of the way: 0 or lin
    if (quadPolyList[i].deg() <= 1) {
        joinPoly(quadPolyList[i], i);  // use a side-effect here: will copy the polynomial to linStore (except: it's zero)
        quadPolyList[i].clear();  // poly is gone -> make space
        emptyRows.push(i);  // also register it as empty
        continue; // next poly - nothing to do here
    }
    // third: re-add the "new" polynomial
    joinPoly(quadPolyList[i], i);
    // check if we are done here
  }

  assert(checkConsistency());  // make sure everything is still fine
  monSet.clear();
  polySet.clear();
  return retNum;
}

/////////////////////////////////////////////////////////////////////
set<uint32> cvEraseVarSet;
bool quadStore_t::commonVars(size_t startVal, size_t endVal, set<uint32> &resSet,
                  monVecQuad_t workingList, unordered_set<monom_t, hash_monom_t> &monSet) {
  if (startVal > endVal) return true;
  // test middle 
  size_t midPos = (endVal+startVal) / 2;

  //linear monomials allow every variable that leads to a monomial in monSet
  monom_t curMon = workingList[midPos];
  cvEraseVarSet.clear();
  // loop over 1..2 elements
  for (auto varItr = resSet.begin(); varItr!=resSet.end(); ++varItr) {
    uint32 var = *varItr;
    // linear or quadratic monomial
    if (monDeg(curMon) <= 1) {   
      monom_t newMon; newMon[0] = var; newMon[1] = curMon[1]; normalizeMonom(newMon);
      if (monSet.find(newMon) == monSet.end()) cvEraseVarSet.insert(var); 
    } else 
      // quadratic case
      if (!(var == curMon[0] || var == curMon[1])) cvEraseVarSet.insert(var); 
  }
  // delete all elements we need to delete
  for (auto itr = cvEraseVarSet.begin(); itr != cvEraseVarSet.end(); ++itr) resSet.erase(*itr);
  // return if we cannot find a set with more than 0 elements
  if (resSet.size() == 0) return false;

  // recursive call
  if (!commonVars(midPos+1, endVal, resSet, workingList, monSet)) return false;
  return commonVars(startVal, midPos-1, resSet, workingList, monSet);
}

// SL for the System 
// its like XL but preserves the monomial structure of the system
size_t quadStore_t::sl() {
  //get the number of equations in the original system
  size_t sysNum = replaceLin.rank() + quadPolyList.size() - emptyRows.size();
  //setting up sets of monomials and vars
  unordered_set< monom_t, hash_monom_t > monSet;
  uint32 otherVar = 0;
  
  for (auto itr = quadPolyByMon.begin(); itr != quadPolyByMon.end(); ++itr)
    monSet.insert((*itr).first);
  for (auto itr = replaceLin.polyByMon.begin(); itr != replaceLin.polyByMon.end(); ++itr){
    monom_t curMon;
    curMon[0] = otherVar;
    curMon[1] = (*itr).first;
    monSet.insert(curMon);
  }
  
  //iterate over all polynomials 
  size_t maxPoly = quadPolyList.size();
  for (size_t curPos = 0; curPos < maxPoly; curPos++)  {
    polyGF2 curPoly = quadPolyList[curPos];
    if (curPoly.isZero() or (curPoly.len() == 1)) continue;
    monVecQuad_t curPolyMon = curPoly.mons;
    
    //looking for suitable equations
    set< uint32 > varSet; varSet.clear(); 
    varSet.insert(curPoly.mons[0][0]); varSet.insert(curPoly.mons[0][1]);
    if (!commonVars(1, curPoly.len()-1, varSet, curPoly.mons, monSet)) continue;

    //calculate the new polynomials and insert them into the system of equations
    monVecQuad_t newMonVec;
    polyGF2 newPoly;
    for (auto varItr = varSet.begin(); varItr!=varSet.end(); ++varItr) {
      newMonVec.clear();
      newPoly.clear();
      uint32 var = *varItr;
      vector<monom_t> curMonSet = curPoly.mons;
      //iterate over all monomials
      for (auto monIter = curMonSet.begin(); monIter != curMonSet.end(); ++monIter) {
        monom_t curMon = *monIter;
	if (monDeg(curMon) == 0) 
	  // deg 0 case
	  curMon[0] = var;
	else if (monDeg(curMon) == 1) {   
	  // deg 1 case
	  assert(curMon[0] == 0);  // have a 1-factor here	  
	  curMon[0] = var;
	} else {
	  // deg 2 case
	  assert(monDeg(curMon) == 2);
	  assert((var == curMon[0]) or (var == curMon[1]));
	}
	// add new monomial to the overall polynomial
	normalizeMonom(curMon);
	newMonVec.push_back(curMon);
      }
      normalizeMons(newMonVec);
      newPoly.mons.swap(newMonVec);
      joinPoly(newPoly);
    } // end of while loop
  }
  //if something has changed start the solver again
  size_t newSysNum = replaceLin.rank() + quadPolyList.size() - emptyRows.size();
  return newSysNum - sysNum; // return the number of newly created polynomials
}

/////////////////////////////////////////////////////////////////////
bool quadStore_t::halfCommonVars(size_t startVal, size_t endVal, set<uint32> &resSet,
                  monVecQuad_t &workingList) {
  if (startVal > endVal) return true;
  // test middle 
  size_t midPos = (endVal+startVal) / 2;

  //linear monomials allow every variable that leads to a monomial in monSet
  monom_t curMon = workingList[midPos];
  cvEraseVarSet.clear();
  // loop over 1..2 elements
  for (auto varItr = resSet.begin(); varItr!=resSet.end(); ++varItr) {
    uint32 var = *varItr;
    // quadratic case
    if ((monDeg(curMon) == 2) and (var != curMon[0]) and (var != curMon[1])) cvEraseVarSet.insert(var); 
  }
  // delete all elements we need to delete
  for (auto itr = cvEraseVarSet.begin(); itr != cvEraseVarSet.end(); ++itr) resSet.erase(*itr);
  // return if we cannot find a set with more than 0 elements
  if (resSet.size() == 0) return false;

  // recursive call
  if (!halfCommonVars(midPos+1, endVal, resSet, workingList)) return false;
  return halfCommonVars(startVal, midPos-1, resSet, workingList);
}

// half SL for the System 
// its like XL but preserves degree two equations
size_t quadStore_t::hsl(float density, float densityMax) {
  //get the number of equations in the original system
  size_t sysNum = replaceLin.rank() + quadPolyList.size() - emptyRows.size();
  
  //iterate over all polynomials 
  size_t maxPoly = quadPolyList.size();
  for (size_t curPos = 0; curPos < maxPoly; curPos++)  {
    polyGF2 curPoly = quadPolyList[curPos];
    if (curPoly.isZero() or (curPoly.len() == 1)) continue;
    
    //looking for suitable equations
    set< uint32 > varSet; varSet.clear(); 
    varSet.insert(curPoly.mons[0][0]); varSet.insert(curPoly.mons[0][1]);
    if (!halfCommonVars(1, curPoly.len()-1, varSet, curPoly.mons)) continue;
    //calculate the new polynomials and insert them into the system of equations
    monVecQuad_t newMonVec;
    polyGF2 newPoly;
    for (auto varItr = varSet.begin(); varItr!=varSet.end(); ++varItr) {
      newMonVec.clear();
      newPoly.clear();
      uint32 var = *varItr;
      vector<monom_t> curMonSet = curPoly.mons;
      //iterate over all monomials
      for (auto monIter = curMonSet.begin(); monIter != curMonSet.end(); ++monIter) {
        monom_t curMon = *monIter;
	if (monDeg(curMon) == 0) 
	  // deg 0 case
	  curMon[0] = var;
	else if (monDeg(curMon) == 1) {   
	  // deg 1 case
	  assert(curMon[0] == 0);  // have a 1-factor here	  
	  curMon[0] = var;
	} else {
	  // deg 2 case
	  assert(monDeg(curMon) == 2);
	  assert((var == curMon[0]) or (var == curMon[1]));
	}
	// add new monomial to the overall polynomial
	normalizeMonom(curMon);
	newMonVec.push_back(curMon);
      }
      normalizeMons(newMonVec);
      newPoly.mons.swap(newMonVec);
      joinPoly(newPoly);
    } // end of while loop
  }
  //echelonize the system to get rid of linear dependent equations 
  if (density < densityMax) echelonizeSparse();
  else echelonizeDense();
  //if something has changed start the solver again
  size_t newSysNum = replaceLin.rank() + quadPolyList.size() - emptyRows.size();
  return newSysNum - sysNum; // return the number of newly created polynomials
}

// counts the total number of monomials in the store
// workload: linear in the number of polynomials
size_t quadStore_t::monNum() {
  size_t res = 0;
  for(auto itr = quadPolyList.begin(); itr != quadPolyList.end(); ++itr) 
    res += (*itr).len(); // get the number of monomials for this polynomial
  return res;
}


// verifies if the data structures are fully consistent with each other
bool quadStore_t::checkConsistency() {

  // count the number of non-zero polynomials in polyList
  size_t cntRes = 0;
  for (size_t i = 0; i < quadPolyList.size(); i++)  
    if (not(quadPolyList[i].isZero())) cntRes++;
  if (cntRes != polyCnt()) {
    cout << "Number of polynomials differs: cnt:cmp is " << cntRes << ":" << polyCnt() << endl;
    return false;
  }

  // check that all keys in "quadPolyByMon" are in fact degree 2
  for (auto pairIter = quadPolyByMon.begin(); pairIter != quadPolyByMon.end(); ++pairIter) {
    if (monDeg(pairIter->first) != 2) {
      cout << "got wrong monomial in quadPolyByMon: " << pairIter->first[1] << ":" << pairIter->first[0] << endl;
      return false;
    } 
  }

  // loop over all polynomials and check that they are registered in all sets
  for (size_t i = 0; i < quadPolyList.size(); i++) {
    if (quadPolyList[i].isZero()) continue; // nothing to do here
    if (not(quadPolyList[i].isSorted())) {
      cout << "polynomial " << i << " is not sorted" << endl;
      quadPolyList[i].printIt();
      return false;
    }
    // loop over all monomials
    for (auto monIter = quadPolyList[i].mons.begin(); monIter != quadPolyList[i].mons.end(); ++monIter) {
      // check in polyByMon (deg 2 monomials)
      if (monDeg(*monIter) == 2) {
        auto monSet = quadPolyByMon[*monIter];
        if (monSet.find(i) == monSet.end()) {
          cout << "Did not find entry " << i << " for monomial " << (*monIter)[1] << ":" << (*monIter)[0] << endl;
          return false;
        }
      }
      // check in polyByVar
      if ((*monIter)[1] == 0) continue; // nothing to do
      auto varSet = quadPolyByVar[(*monIter)[1]];
      if (varSet.find(i) == varSet.end()) {
        cout << "Did not find entry " << i << " for variable 1 in " << (*monIter)[1] << ":" << (*monIter)[0] << endl;
        return false;
      }
      if ((*monIter)[0] == 0) continue; // nothing to do
      varSet = quadPolyByVar[(*monIter)[0]];
      if (varSet.find(i) == varSet.end()) {
        cout << "Did not find entry " << i << " for variable 0 in " << (*monIter)[1] << ":" << (*monIter)[0] << endl;
        return false;
      }
    }
  }

  // check if set2poly is correct
  // to do

  // check that all empty lines are in fact empty (and vice versa)
  freePos_t restack; // needed to access emptyRows
  unordered_set<int> reSet(emptyRows.size());
  while (not(emptyRows.empty())) {
    int tmpVar = emptyRows.top(); 
    restack.push(tmpVar);
    reSet.insert(tmpVar);
    emptyRows.pop();
  }
  while (not(restack.empty())) {
    int curLine = restack.top();
    restack.pop();
    emptyRows.push(curLine);
    // check if the line is actually empty
    if (!quadPolyList[curLine].isZero()) {
      cout << "non-empty row " << curLine << " with poly ";
      quadPolyList[curLine].printIt();
      return false;
    }
  }
  // loop quadPolyList to check if an empty row is actually in emptyRows (aka reSet)
  for (size_t i = 0; i < quadPolyList.size(); i++) {
    if ((quadPolyList[i].isZero()) and !(reSet.find(i) != reSet.end())) {
      cout << "empty row " << i << " that is not in emptyRows" << endl;
      return false;
    }
  }

  return true; // keep the compiler happy
}


// compute some statistics about the shortest linear chains that could be replaced
void quadStore_t::replaceStat(uint32 &cntVal, uint32 &minVal, uint32 &maxVal, uint32 &medVal, float &avgVal) {
  // collect a list of all variables that may be replaced
  vector<uint32> lenVec;
  lenVec.reserve(quadPolyByVar.size());
  for (auto pairIter = quadPolyByVar.begin(); pairIter != quadPolyByVar.end(); ++pairIter) {
    uint32 varNum = pairIter->first;
    // check if we know something about this variable in linStore
    int varVal = replaceLin.giveValue(varNum); 
    if (varVal == giveValue_unknown) continue;
    // clear value -> length 0
    if ((varVal == giveValue_zero) or (varVal == giveValue_one)) {
      lenVec.push_back(0);
      continue;
    }
    // only one case left - giveValue_linear
    assert(varVal == giveValue_linear);
    lenVec.push_back(replaceLin.linPolyMat.find(varNum)->second.len()); 
  }
 
  // get a silly case out of the way
  if (lenVec.size() == 0) {
    cntVal = minVal = maxVal = medVal = 0;
    avgVal = 0;
    return;
  }  

  // extract all values from this list
  sort(lenVec.begin(), lenVec.end());
  cntVal = lenVec.size();
  minVal = lenVec[0]; maxVal = lenVec[lenVec.size()-1];
  medVal = lenVec[ lenVec.size() / 2 ]; // near-median ;-)
  size_t sumVal = 0;
  for (auto lenIter = lenVec.begin(); lenIter != lenVec.end(); ++lenIter) 
    sumVal += *lenIter;
  avgVal = sumVal;
  avgVal /= lenVec.size();
}

// outputs the full list of polynomials
void quadStore_t::printIt(bool withNum = false) {
  cout << "number quadratic polynomials: " << polyCnt() << endl;
  for(size_t pos = 0; pos < quadPolyList.size(); pos++) {
    if (quadPolyList[pos].isZero()) continue;
    if (withNum) cout << pos << ": ";
    quadPolyList[pos].printIt(); 
  }
}


//////////////////////////////////////////////////////////////////////////////////////////////////////
quadStore_t quadStore;

//////////////////////////////////////////////////////////////////////////////////////////////////////
int main() {
  // some initial testing
  /*
  const int MAXITEMS = 4;
  string allVars[MAXITEMS] = {
    "x0", "x1", "x2", "x3",
  };

  for (int i = 0; i < MAXITEMS; i++) {
    cout << "add vars: " << allVars[i] << "=" << fullVarDict.giveNum(allVars[i]) << endl; 
  }
  
  string inStr("x0*x1 + x1*x2 + x0 + x3");
  polyGF2 inPoly(inStr);
  cout << "add poly " << inStr << " / became: "; inPoly.printIt();

  return 0;
  */
  /*
  const int MAXITEMS = 17;
  string allPolys[MAXITEMS] = {
    "x000+x000", "x001+x001", "x002+x002", "x003+x003", "x004+x004",
    "x000*x001", "x001*x002 + x001*x003", "x000*x001 + x000*x002", "x000*x001 + x000*x002 + x003", "x001*x002 + x002*x003 + x000*x004", 
    "x000*x002 + x001*x003 + x001*x004 + x000", "x000*x002 + x000*x003 + x003*x004 + x001 + x002 + 1", "x001*x002 + x001*x003",
    "x001*x002 + x003*x004 + x003", "x001*x002 + x000*x003 + x002*x003 + x000*x004 + x001*x004", "x001*x004 + x003",  
    "x000*x001 + x001*x002 + x003*x004 + x001 + 1"
  };
  */
  /*
  const int MAXITEMS = 11;
  string allPolys[MAXITEMS] = {
    "a+a", "b+b", "c+c", "d+d",
    "a+b+c", "a*b + c*d + d", "a*b + c*a + a", "c*a + c*b +b", "c*d + a*c + a*b + c", "c*a + a*b + a", "a*c + a + b + c"
  };
  */
  /*
  for (int i = 0; i < MAXITEMS; i++) {
    polyGF2 inPoly(allPolys[i]);
    cout << "just add: " << allPolys[i] << endl; 
    cout << "  became: "; inPoly.printIt();
    quadStore.insertPoly(inPoly);
  }

  string curVar;
  curVar = "x000"; cout << "var: " << curVar << " is " << fullVarDict.giveNum(curVar) << endl;
  curVar = "x001"; cout << "var: " << curVar << " is " << fullVarDict.giveNum(curVar) << endl;
  curVar = "x002"; cout << "var: " << curVar << " is " << fullVarDict.giveNum(curVar) << endl;
  curVar = "x003"; cout << "var: " << curVar << " is " << fullVarDict.giveNum(curVar) << endl;
  curVar = "x004"; cout << "var: " << curVar << " is " << fullVarDict.giveNum(curVar) << endl;
  curVar = "x005"; cout << "var: " << curVar << " is " << fullVarDict.giveNum(curVar) << endl;

  for (uint i = 0; i < 3; i++) {

    cout << ">>>result<<<" << endl;

    quadStore.printIt();
    cout << endl;
    replaceLin.printIt();

    //cout << "soft reduce" << endl;
    //quadStore.softReduce(1 << 31, 5*i);

    quadStore.printIt();
    cout << endl;
    replaceLin.printIt();

    cout << "echelonize sparse" << endl;
    quadStore.echelonizeSparse();

    cout << "eval: " << quadStore.fullEval() << endl;
  
    //quadStore.printIt();
    //cout << endl;
    //replaceLin.printIt();

    cout << ">>>result<<<" << endl;

    quadStore.printIt();
    cout << endl;
    replaceLin.printIt();

  }

  return 0;

  cout << endl;
  cout << "echelonize" << endl;
  quadStore.echelonizeDense();

  cout << ">>>result<<<" << endl;

  quadStore.printIt();
  cout << endl;
  replaceLin.printIt();


  return 0;

  quadStore.printIt();
  cout << endl;
  replaceLin.printIt();
  
  cout << "eval: " << quadStore.fullEval() << endl;

  quadStore.printIt();
  cout << endl;
  replaceLin.printIt();

  cout << "eval: " << quadStore.fullEval() << endl;

  quadStore.printIt();
  cout << endl;
  replaceLin.printIt();

  cout << endl;
  // get all numbers
  for (auto itr = fullVarDict.num2var.begin(); itr != fullVarDict.num2var.end(); ++itr) {
    cout << itr->first << ":" << itr->second << "  ";
  }
  cout << endl;

  return 0;
*/

  // check the predefined types; just a precaution
  if (sizeof(uint32) != 4) {
    cout << "Error: Underlying datatype has wrong size, namely " << sizeof(uint32) << " byte" << endl;
    return 1;
  }
  assert(sizeof(uint32) == 4); // in case we have "assert"

  // initialize variables
  string connectionID, command;
  stringstream convert;
  string curLine;
  clock_t timeSL, timeDense, timeSparse, timeEval;
  timeSL = timeDense = timeSparse = timeEval = 0;

  // protocol: 
  // "connectionID command parameters \n", e.g.
  // C> 123 checkConnection
  // <S 123 check connection O.K.
  // C> 345 end
  // <S 345 terminated
  while (true) {
    bool knownCmd = false;
    // number to identify this connection
    cin >> connectionID;
    // get the command
    cin >> command;

    // debug only
    //cout << "got command: " << command << endl;

    if (command == "end") {
      knownCmd = true;
      cout << connectionID << " terminated" << endl;
      break; // quit loop
    }

    if (command == "checkConnection") {
      knownCmd = true;
      // just mirror the command
      cout << connectionID << " check connection O.K." << endl;
    }

    if (command == "addEqs") {
      knownCmd = true;
      // read in all polynomials
      //cout << "got: " << curLine << endl;
      bool keepReading = true;
      while (keepReading) { 
        getline(cin, curLine);
        // check if there is any data at all
        bool noData = true; size_t strPos = 0;
        while (noData and (strPos < curLine.size())) {
          // check if we have variables here
          if (((curLine[strPos] >= 'a') and (curLine[strPos] <= 'z')) or
              ((curLine[strPos] >= 'A') and (curLine[strPos] <= 'Z')) or
              ((curLine[strPos] >= '0') and (curLine[strPos] <= '9')) or
              (curLine[strPos] == '_')) noData = false;
          // check for the end character
          if (curLine[strPos] == '=') { 
            // we are done here ('=')
            noData = true; keepReading = false;
            strPos = curLine.size();
            break;
          }
          strPos++; 
        }
        if (noData) continue;
        polyGF2 inPoly(curLine);
        quadStore.insertPoly(inPoly);
      } // end while (curLine != "-1");
      // quit O.K.
      cout << connectionID << " add eqs O.K." << endl;
    }

    if (command == "getVars") {
      knownCmd = true;
      // follow the protocol
      cout << connectionID << " ";
      getline(cin, curLine);
      // get all variables and reflect their value (if any) to cout
      bool isFirst = true;
      if (not(replaceLin.isConsistent)) {
        cout << "!not consistent!=1";
        isFirst = false;
      }
      string curVar = "";
      for (size_t i = 0; i <= curLine.size(); i++) {
        char curChar;
        if (i < curLine.size()) curChar = curLine[i]; else curChar = ' ';
        if ((curChar != ' ') and (i < curLine.size())) { curVar += curChar; continue; }
        if (curVar.size() == 0) continue;
        // check the value of this variable
        int varVal = replaceLin.giveValue(fullVarDict.giveNum(curVar));
        if ((varVal == 0) or (varVal == 1)) {
          if (not(isFirst)) cout << ",";
          cout << curVar << "=" << varVal;
          isFirst = false;
        }
        curVar = "";
      }
      cout << endl;
    }

    if (command == "statChains") {
      knownCmd = true;
      // checks how long the quadratic & linear equations are
      // start with the quadratic equations
      vector<int> sizeVec;
      sizeVec.resize(100); // just add some zeros
      for (auto polyIt = quadStore.quadPolyList.begin(); polyIt != quadStore.quadPolyList.end(); ++polyIt) {
        size_t len = polyIt->len();
        if (len == 0) continue; // skip empty lines
        while (sizeVec.size() <= len) sizeVec.resize(2*len);
        sizeVec[len]++;
      }
      // get the result to the caller
      cout << connectionID << "+quad " ;
      bool first = true;
      for (size_t i = 0; i < sizeVec.size(); i++) {
        if (sizeVec[i] == 0) continue;
        if (not(first)) cout << ",";
        cout << i << ":" << sizeVec[i];
        first = false;
      } 
      cout << endl;

      // get information about the linear equations
      sizeVec.clear();
      for (auto polyIt = replaceLin.linPolyMat.begin(); polyIt != replaceLin.linPolyMat.end(); ++polyIt) {
        size_t len = polyIt->second.len();
        if (len == 0) continue; // skip empty lines
        if (sizeVec.size() <= len) sizeVec.resize(2*len);
        sizeVec[len]++;
      }
      // get the result to the caller
      cout << connectionID << " lin " ;
      first = true;
      for (size_t i = 0; i < sizeVec.size(); i++) {
        if (sizeVec[i] == 0) continue;
        if (not(first)) cout << ",";
        cout << i << ":" << sizeVec[i];
        first = false;
      } 
      cout << endl;
    }


    if (command == "getLinSys") {
      knownCmd = true;
      if (not(replaceLin.isConsistent)) cout << connectionID << "+" << "!not consistent!" << endl;
      // output the full system
      linPolyMat_t::const_iterator itr;
      for(itr = replaceLin.linPolyMat.begin(); itr != replaceLin.linPolyMat.end(); ++itr) {
        linPoly ourPoly = (*itr).second;
        cout << connectionID << "+";
        ourPoly.printIt(); 
      } 
      cout << connectionID << " 0" << endl; // last polynomial
    }

    if (command == "getQuadSys") {
      knownCmd = true;
      if (not(replaceLin.isConsistent)) cout << connectionID << "+" << "!not consistent!" << endl;
      // output the full system
      for(auto itr = quadStore.quadPolyList.begin(); itr != quadStore.quadPolyList.end(); ++itr) {
        if (itr->isZero()) continue;
        cout << connectionID << "+";
        itr->printIt(); 
      } 
      cout << connectionID << " 0" << endl; // last polynomial
    }

    if (command == "fullEval") {
      knownCmd = true;
      size_t maxMonNum;
      cin >> maxMonNum;
      size_t outMax;
      cin >> outMax;
      size_t inMax;
      cin >> inMax;

      // return the number of equations changed
      clock_t startTime = clock();
      cout << connectionID << " fullEval changed " << quadStore.fullEval(maxMonNum, outMax, inMax) << endl;
      timeEval += clock()-startTime;
    }

    if (command == "fullSparseEval") {
      knownCmd = true;
      // return the number of equations changed
      clock_t startTime = clock();
      cout << connectionID << " fullSEval changed " << quadStore.fullSparseEval() << endl;
      timeEval += clock()-startTime;
    }

                    
    if (command == "echelonizeSparse") {
      knownCmd = true;
      // get the maximal size
      size_t maxSize;
      cin >> maxSize;
      clock_t startTime = clock();
      quadStore.echelonizeSparse(maxSize);
      timeSparse += clock()-startTime;
      // return the number of equations changed
      cout << connectionID << " done echelonizeSparse" << endl;
    }

    if (command == "echelonizeDense") {
      knownCmd = true;
      clock_t startTime = clock();
      quadStore.echelonizeDense();
      timeDense += clock()-startTime;
      // return the number of equations changed
      cout << connectionID << " done echelonizeDense" << endl;
    }

    if (command == "sl") {
      knownCmd = true;
      // return the number of equations changed
      clock_t startTime = clock();
      cout << connectionID << " sl changed " << quadStore.sl() << endl;
      timeSL += clock()-startTime;
    }

    if (command == "hsl") {
      knownCmd = true;
      float curSize;
      float maxSize;
      cin >> curSize;
      cin >> maxSize;
      // return the number of equations changed
      clock_t startTime = clock();
      cout << connectionID << " hsl changed " << quadStore.hsl(curSize,maxSize) << endl;
      timeSL += clock()-startTime;
    }

    if (command == "softReduce") {
      knownCmd = true;
      size_t maxMonNum;
      cin >> maxMonNum;
      size_t maxPolyExpand;
      cin >> maxPolyExpand;
      // return the number of equations changed
      quadStore.softReduce(maxMonNum, maxPolyExpand);
      cout << connectionID << " done softReduce" << endl;
    }

    if (command == "stats") {
      knownCmd = true;
      unsigned int level;
      cin >> level;
      // return all relevant stats we have for this level
      if (level >= 1) {
        // all data together
        cout << connectionID << "+allVars " << fullVarDict.varsNum() << endl;
        cout << connectionID << "+allRows " << replaceLin.rank()+quadStore.polyCnt() << endl;

        // linear data
        cout << connectionID << "+linIsConsistent " << (int)replaceLin.isConsistent << endl;
        cout << connectionID << "+linRows " << replaceLin.rank() << endl;
        cout << connectionID << "+linCols " << replaceLin.polyByMon.size() << endl; // also the number of variables in use

        // quadratic data
        cout << connectionID << "+quadVarsInUse " << quadStore.quadPolyByVar.size() << endl;
        cout << connectionID << "+quadRows " << quadStore.polyCnt() << endl;
        cout << connectionID << "+quadCols " << quadStore.quadPolyByMon.size() << endl;

        // timing infos
        cout << connectionID << "+time.sl " << timeSL << endl;
        cout << connectionID << "+time.dense " << timeDense << endl;
        cout << connectionID << "+time.sparse " << timeSparse << endl;
        cout << connectionID << "+time.eval " << timeEval << endl;
        cout << connectionID << "+time.cps " << CLOCKS_PER_SEC << endl;
      }
      if (level >= 2) {
        cout << connectionID << "+linNonZeroEntries " << replaceLin.monNum() << endl;
        cout << connectionID << "+quadNonZeroEntries " << quadStore.monNum() << endl;
        uint32 cntVal, minVal,maxVal,medVal;
        float avgVal;
        quadStore.replaceStat(cntVal, minVal,maxVal,medVal,avgVal);
        cout << connectionID << "+allReplaceChains.cnt " << cntVal << endl;
        cout << connectionID << "+allReplaceChains.min " << minVal << endl;
        cout << connectionID << "+allReplaceChains.max " << maxVal << endl;
        cout << connectionID << "+allReplaceChains.med " << medVal << endl;
        cout << connectionID << "+allReplaceChains.avg " << avgVal << endl;
      }
      cout << connectionID << " done stats" << endl;
    }

    if (not(knownCmd)) {
      cout << connectionID << " error: unknown command " << command << endl;
    }
  } // end while loop
  return 0; // we are done
}

