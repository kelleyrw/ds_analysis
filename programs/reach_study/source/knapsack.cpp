#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include "localsolver.h"
#include <boost/numpy.hpp>
#include <boost/scoped_array.hpp>

using namespace localsolver;
using namespace std;
namespace bp = boost::python;
namespace np = boost::numpy;

class Knapsack {
    public:

        // constructor
        Knapsack
        (
            np::ndarray const & np_weights,
            np::ndarray const & np_values,
            int const weight_limit
        )
            : nbItems(np_weights.shape(0))
            , weights(np_weights) 
            , values(np_values) 
            , knapsackBound(weight_limit)
            , solution(np::empty(np_weights.get_nd(), np_weights.get_shape(), np_weights.get_dtype()))
        {
        }

        /* Number of items. */
        int nbItems;

        /* Items properties. */
        np::ndarray weights;
        np::ndarray values;

        /* Knapsack bound */
        lsint knapsackBound;

        /* LocalSolver. */
        LocalSolver localsolver;

        /* Decision variables. */
        std::vector<localsolver::LSExpression> x;

        /* Solution (items in the knapsack). */
        //std::vector<int> solution;
        np::ndarray solution;


        // envoke the solver
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
                    LSExpression itemWeight = model.createExpression(O_Prod, x[i], static_cast<lsint>( bp::extract<int>(weights[i])));
                    weightSum.addOperand(itemWeight);
                }    
                LSExpression weightConstraint = model.createExpression(O_Leq, weightSum, knapsackBound);
                model.addConstraint(weightConstraint);

                // maximize value
                LSExpression valueSum = model.createExpression(O_Sum);
                for (int i = 0; i < nbItems; i++) {
                    LSExpression itemValue = model.createExpression(O_Prod, x[i], static_cast<lsint>( bp::extract<int>(values[i])));
                    valueSum.addOperand(itemValue);
                }
                model.addObjective(valueSum, OD_Maximize);
                model.close();

                LSPhase phase = localsolver.createPhase();
                phase.setTimeLimit(limit);
                localsolver.solve();

                //solution.clear();
                for (int i = 0; i < nbItems; ++i)
                {
                    solution[i] = static_cast<int>(x[i].getValue());
                }

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

            for (unsigned int i = 0; i < nbItems; ++i)
                outfile << bp::extract<int>(solution[i]) << " ";
            outfile << endl;
            outfile.close();
        }

        void printSolution() 
        {
            for (unsigned int i = 0; i < nbItems; ++i)
            {
                std::cout << bp::extract<int>(solution[i]) << " ";
            }
            std::cout << std::endl;
        }

};

np::ndarray knapsack_solve
(
    np::ndarray const & weights,
    np::ndarray const & values,
    int const weight_limit, 
    int const time_limit 
) 
{
    Knapsack model(weights, values, weight_limit);
    model.solve(time_limit);
    //model.writeSolution("output/sol.txt");
    return model.solution;
}

BOOST_PYTHON_MODULE(knapsack) {
    np::initialize();  // have to put this in any module that uses Boost.NumPy
    bp::def("solve", knapsack_solve);
}
