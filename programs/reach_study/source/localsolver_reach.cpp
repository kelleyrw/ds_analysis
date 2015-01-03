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

class OptimizeTRP {
    public:

        // constructor
        OptimizeTRP
        (
            np::ndarray const & np_grps,
            np::ndarray const & np_trps,
            double const bp_grp_limit
        )
            : num_spots(np_grps.shape(0))
            , grps(np_grps) 
            , trps(np_trps) 
            , grp_limit(bp_grp_limit)
            , solution(np::empty(np_grps.shape(0), 
                                 np_grps.get_shape(), 
                                 np::dtype::get_builtin<int>()))
        {
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
        void solve(int time_limit) {
            try {
                LSModel model = localsolver.getModel();

                // decision variables u[i]
                u.resize(num_spots);
                for (int i = 0; i < num_spots; i++) {
                    u[i] = model.createExpression(O_Int, lsint(0), lsint(10));
                    stringstream s;
                    s << "u[" << i << "]";
                    u[i].setName(s.str());
                }

                // GRP constraint
                LSExpression grp_sum = model.createExpression(O_Sum);
                for (int i = 0; i < num_spots; i++) {
                    LSExpression item_grp = model.createExpression(O_Prod, u[i], static_cast<lsdouble>( bp::extract<double>(grps[i])));
                    grp_sum.addOperand(item_grp);
                }    
                LSExpression grp_constraint = model.createExpression(O_Leq, grp_sum, grp_limit);
                model.addConstraint(grp_constraint);

                // maximize trp
                LSExpression trp_sum = model.createExpression(O_Sum);
                for (int i = 0; i < num_spots; i++) {
                    LSExpression item_trp = model.createExpression(O_Prod, u[i], static_cast<lsdouble>( bp::extract<double>(trps[i])));
                    trp_sum.addOperand(item_trp);
                }
                model.addObjective(trp_sum, OD_Maximize);
                model.close();

                LSPhase phase = localsolver.createPhase();
                phase.setTimeLimit(time_limit);
                localsolver.solve();

                //solution.clear();
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
    int const time_limit 
) 
{
    OptimizeTRP model(grps, trps, grp_limit);
    model.solve(time_limit);
    model.printSolution();
    return model.solution;
}

BOOST_PYTHON_MODULE(localsolver_reach) {
    np::initialize();  // have to put this in any module that uses Boost.NumPy
    bp::def("optimize", optimize);
}
