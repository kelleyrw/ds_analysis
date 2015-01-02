#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include "localsolver.h"

using namespace localsolver;
using namespace std;

class Reach {
public:
  /* Number of items. */
  int nbItems;

  /* Items properties. */
  std::vector<lsint> weights;
  std::vector<lsint> values;

  /* Knapsack bound */
  lsint knapsackBound;
  
  /* LocalSolver. */
  LocalSolver localsolver;

  /* Decision variables. */
  std::vector<localsolver::LSExpression> x;

  /* Solution (items in the knapsack). */
  std::vector<int> solution;

  void readInstance(const string& fileName) {
    ifstream infile(fileName.c_str());
    if (!infile.is_open()) {
      cerr << "File " << fileName << " cannot be opened." << endl;
      exit(1);
    }
    infile >> nbItems;

    weights.resize(nbItems);
    for (int i = 0; i < nbItems; i++)
      infile >> weights[i];

    values.resize(nbItems);
    for (int i = 0; i < nbItems; i++)
      infile >> values[i];
  
    infile >> knapsackBound;
  
    infile.close();
    cout << "Input read" << endl;
  }

  void solve(int limit) {
    try {
      LSModel model = localsolver.getModel();

      // decision variables x[i] 
      x.resize(nbItems);
      for (int i = 0; i < nbItems; i++) {
        x[i] = model.createExpression(O_Bool);
        stringstream s;
        s << "x[" << i << "]";
        x[i].setName(s.str());
      }

      // weight constraint
      LSExpression weightSum = model.createExpression(O_Sum);
      for (int i = 0; i < nbItems; i++) {
        LSExpression itemWeight = model.createExpression(O_Prod, x[i], weights[i]);
        weightSum.addOperand(itemWeight);
      }    
      LSExpression weightConstraint = model.createExpression(O_Leq, weightSum, knapsackBound);
      model.addConstraint(weightConstraint);
    
      // maximize value
      LSExpression valueSum = model.createExpression(O_Sum);
      for (int i = 0; i < nbItems; i++) {
        LSExpression itemValue = model.createExpression(O_Prod, x[i], values[i]);
        valueSum.addOperand(itemValue);
      }
      model.addObjective(valueSum, OD_Maximize);
      model.close();

      LSPhase phase = localsolver.createPhase();
      phase.setTimeLimit(limit);
      localsolver.solve();

      solution.clear();
      for (int i = 0; i < nbItems; ++i)
        if (x[i].getValue() == 1) 
          solution.push_back(i);
		  
    } catch (LSException *e) {
      cout << "LSException:" << e->getMessage() << std::endl;
      exit(1);
    }
  }

  void writeSolution(const string& fileName) {
    ofstream outfile(fileName.c_str());
    if (!outfile.is_open()) {
      cerr << "File " << fileName << " cannot be opened." << endl;
      exit(1);
    }

    for (unsigned int i = 0; i < solution.size(); ++i)
      outfile << solution[i] << " ";
    outfile << endl;
    outfile.close();
  }
};

int main(int argc, char** argv) {
  if (argc != 2) {
    cout << "Usage: knapsack inputFile" << endl;
    exit(1);
  }
  char *instanceFile = argv[1];  

  Knapsack model;
  model.readInstance(instanceFile);
  model.solve(10);
  model.writeSolution("output/sol.txt");

  return 0;
}

