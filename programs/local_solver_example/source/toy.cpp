/********** toy.cpp **********/
#include <iostream>
#include "localsolver.h"

using namespace localsolver;

int main()
try
{
    lsint weights[] = {10, 60, 30, 40, 30, 20, 20, 2};
    lsint values[] = {1, 10, 15, 40, 60, 90, 100, 15};
    lsint knapsackBound = 102;

    LocalSolver localsolver;
    LSModel model = localsolver.getModel();

    // 0-1 decisions
    LSExpression x[8];
    for (int i = 0; i < 8; i++) 
        x[i] = model.createExpression(O_Bool);

    // knapsackWeight <- 10*x0 + 60*x1 + 30*x2 + 40*x3 + 30*x4 + 20*x5 + 20*x6 + 2*x7;
    LSExpression knapsackWeight = model.createExpression(O_Sum);
    for (int i = 0; i < 8; i++) 
        knapsackWeight.addOperand(model.createExpression(O_Prod, weights[i], x[i]));

    // knapsackWeight <= knapsackBound;
    model.addConstraint(model.createExpression(O_Leq, knapsackWeight, knapsackBound));

    // knapsackValue <- 1*x0 + 10*x1 + 15*x2 + 40*x3 + 60*x4 + 90*x5 + 100*x6 + 15*x7;
    LSExpression knapsackValue = model.createExpression(O_Sum);
    for (int i = 0; i < 8; i++) 
        knapsackValue.addOperand(model.createExpression(O_Prod, values[i], x[i]));

    // maximize knapsackValue;
    model.addObjective(knapsackValue, OD_Maximize);

    // close model, then solve
    model.close();
    LSPhase phase = localsolver.createPhase();
    phase.setTimeLimit(1);
    localsolver.solve();

    return 0;

} 
catch (const LSException& e)
{
    std::cerr << "LSException occured:" << e.getMessage() << std::endl;
    return 1;
}
