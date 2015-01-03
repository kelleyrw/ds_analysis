#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include "localsolver.h"
#include "NumPyArrayData.hpp" 
#include <boost/numpy.hpp>
#include <boost/scoped_array.hpp>

using namespace localsolver;
using namespace std;
namespace bp = boost::python;
namespace np = boost::numpy;

class OptimizeTRP {
    public:

        // constructor
        OptimizeTRP
        (
            np::ndarray const & np_grps,
            np::ndarray const & np_trps,
            double const bp_grp_limit,
            int optimization_option
        )
            : num_spots(np_grps.shape(1))
            , grps(np_grps) 
            , trps(np_trps) 
            , grp_limit(bp_grp_limit)
            , solution(np::zeros(bp::make_tuple(num_spots), 
                                 np::dtype::get_builtin<int>()))
        {
            std::cout << "instantiating the class"  << std::endl;
        }

        /* Number of spots. */
        int num_spots;

        /* Total and Target impressions*/
        np::ndarray grps;
        np::ndarray trps;

        /* GRP bound */
        lsdouble grp_limit;

        /* LocalSolver. */
        LocalSolver localsolver;

        /* Decision variables, i.e. units awarded */
        std::vector<localsolver::LSExpression> u;

        /* Solution (items in the knapsack). */
        //std::vector<int> solution;
        np::ndarray solution;

        // envoke the solver
        void solve(int time_limit, int annealing_level) {
            try {
                LSModel model = localsolver.getModel();

                // ndarray accessor wrappers
                NumPyArrayData<double> grps_data(grps);
                NumPyArrayData<double> trps_data(trps);
                //NumPyArrayData<int>    solution_data(solution);

                // decision variables u[i]
                std::cout << "number of spots: " << num_spots << std::endl;
                u.resize(num_spots);
                for (int i = 0; i < num_spots; i++) {
                    u[i] = model.createExpression(O_Int, lsint(0), lsint(10));
                    stringstream s;
                    s << "u[" << i << "]";
                    u[i].setName(s.str());
                }

                // GRP constraint
                std::cout << "adding constraints" << std::endl;
                LSExpression grp_sum = model.createExpression(O_Sum);
                for (int i = 0; i < num_spots; i++) {
                    LSExpression item_grp = model.createExpression(O_Prod, u[i], lsdouble(grps_data(0,i)));
                    grp_sum.addOperand(item_grp);
                }    
                LSExpression grp_constraint = model.createExpression(O_Leq, grp_sum, grp_limit);
                model.addConstraint(grp_constraint);

                grp_constraint = model.createExpression(O_Geq, grp_sum, grp_limit-1.0);
                model.addConstraint(grp_constraint);

                // maximize trp
                std::cout << "adding objective function" << std::endl;
                LSExpression trp_sum = model.createExpression(O_Sum);
                for (int i = 0; i < num_spots; i++) {
                    LSExpression item_trp = model.createExpression(O_Prod, u[i], lsdouble(trps_data(0,i)));
                    trp_sum.addOperand(item_trp);
                }
                model.addObjective(trp_sum, OD_Maximize);
                model.close();

                // Local Solver tesstings
                LSPhase phase = localsolver.createPhase();
                phase.setTimeLimit(time_limit);

                localsolver.getParam().setAnnealingLevel(annealing_level);
                
                std::cout << "\ncalling local solver" << std::endl;
                localsolver.solve();

                //solution.clear();
                std::cout << "i\nprinting and returning solution" << std::endl;
                for (int i = 0; i < num_spots; ++i)
                {
                    solution[i] = static_cast<int>(u[i].getValue());
                }

            } catch (LSException *e) {
                cout << "LSException:" << e->getMessage() << std::endl;
                exit(1);
            }
        }

        void printSolution() 
        {
            for (unsigned int i = 0; i < num_spots; ++i)
            {
                std::cout << bp::extract<int>(solution[i]) << " ";
            }
            std::cout << std::endl;
        }

};

np::ndarray optimize
(
    np::ndarray const & grps,
    np::ndarray const & trps, 
    double const grp_limit, 
    int const time_limit = 10,
    int optimization_option = 0,
    int const annealing_level = 0
) 
{
    std::cout << "\nSet up local solver" << std::endl;
    OptimizeTRP model(grps, trps, grp_limit, optimization_option);
    std::cout << "\nsolving" << std::endl;
    model.solve(time_limit, annealing_level);
    return model.solution;
}


// ndarray accessor wrappers
//NumPyArrayData<double> grps_data(grps);
//NumPyArrayData<double> trps_data(trps);
void test_numpy_wrapper( np::ndarray input_array)
{
    int num_rows = input_array.shape(0);
    int num_cols = input_array.shape(1);

    NumPyArrayData<double> wrapped_array(input_array);

    std::cout << "\nM = [\t";
    for (int i = 0; i < num_rows; ++i)
    {
        if( i != 0) 
        {
            std::cout << "\t"; 
        }
        for (int j = 0; j < num_cols; ++j)
        {
            std::cout << " " << wrapped_array(i,j);
        }
        if( i == num_rows - 1) 
        {
            std::cout << " ]\n" <<  std::endl;
        } 
        else
        {
            std::cout <<  std::endl;
        }
    }
}
                
                

BOOST_PYTHON_MODULE(localsolver_reach) {
    np::initialize();  // have to put this in any module that uses Boost.NumPy
    bp::def("optimize", optimize);
    bp::def("test", test_numpy_wrapper);
}
