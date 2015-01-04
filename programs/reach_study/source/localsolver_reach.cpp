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
            double const bp_grp_limit
        )
            : num_spots(np_grps.shape(1))
            , grps(np_grps) 
            , trps(np_trps) 
            , grp_limit(bp_grp_limit)
            , solution(np::zeros(bp::make_tuple(num_spots), 
                                 np::dtype::get_builtin<int>()))
        {
            std::cout << "instantiating an Optimize TRP class"  << std::endl;
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


                // maximize for TRP
                LSExpression trp_sum = model.createExpression(O_Sum);
                for (int i = 0; i < num_spots; i++) {
                    LSExpression item_trp = model.createExpression(O_Prod, u[i], lsdouble(trps_data(0,i)));
                    trp_sum.addOperand(item_trp);
                }
                model.addObjective(trp_sum, OD_Maximize);
                model.close();

                // Local Solver settings
                LSPhase phase = localsolver.createPhase();
                phase.setTimeLimit(time_limit);

                // set the annealing level
                localsolver.getParam().setAnnealingLevel(annealing_level);
                
                std::cout << "\ncalling local solver" << std::endl;
                localsolver.solve();

                //solution.clear();
                std::cout << "\nprinting and returning solution" << std::endl;
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

class OptimizeReach {
    public:

        // constructor
        OptimizeReach
        (
            np::ndarray const & np_log_p,  // log(p_iP), where i is the individual, and P is a program
            np::ndarray const & np_w,      // weight of individual i
            np::ndarray const & np_grps,   // gross rating point for "total" audience
            double const bp_grp_limit      // specified limit on delievered GRPs
        )
            : num_individuals(np_log_p.shape(0))
            , num_spots(np_log_p.shape(1))
            , log_p(np_log_p)
            , w(np_w)
            , grps(np_grps) 
            , grp_limit(bp_grp_limit)
            , solution(np::zeros(bp::make_tuple(num_spots), 
                                 np::dtype::get_builtin<int>()))
        {
            std::cout << "instantiating an Optimize Reach class"  << std::endl;
            std::cout << "number of individuals: " << num_individuals << std::endl;
            std::cout << "number of spots: "       << num_spots << std::endl;
            
            std::cout << "dim of w = (" << w.shape(0) << ", " << w.shape(1) << ")" << std::endl;
        }

        /* Number of individuals and spots. */
        int num_individuals;
        int num_spots;

        /* Total impressions in GRPs*/
        np::ndarray log_p;
        np::ndarray w;
        np::ndarray grps;

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
                NumPyArrayData<double> log_p_data(log_p);
                NumPyArrayData<double> w_data(w);
                NumPyArrayData<double> grps_data(grps);

                // decision variables u[i]
                std::cout << "number of individuals: " << num_individuals << std::endl;
                std::cout << "number of spots: "       << num_spots << std::endl;
                u.resize(num_spots);
                for (int P = 0; P < num_spots; P++) {
                    u[P] = model.createExpression(O_Int, lsint(0), lsint(10));
                    stringstream s;
                    s << "u[" << P << "]";
                    u[P].setName(s.str());
                }

                // GRP constraint
                std::cout << "adding constraints" << std::endl;
                LSExpression grp_sum = model.createExpression(O_Sum);
                for (int P = 0; P < num_spots; P++) {
                    LSExpression item_grp = model.createExpression(O_Prod, u[P], lsdouble(grps_data(0,P)));
                    grp_sum.addOperand(item_grp);
                }    
                LSExpression grp_constraint = model.createExpression(O_Leq, grp_sum, grp_limit);
                model.addConstraint(grp_constraint);

                grp_constraint = model.createExpression(O_Geq, grp_sum, grp_limit-1.0);
                model.addConstraint(grp_constraint);


                // maximize for Reach
                // build the target universe
                double target_universe = 0.0; 

                // buld outer summation over individuals
                LSExpression sum_Individuals = model.createExpression(O_Sum);
                for (int i = 0; i < num_individuals; i++) 
                {
                    // update universe
                    target_universe += w_data(i,0);

                    // for each individual, loop over the programs P
                    // build argument of the exponential
                    LSExpression sum_logPU = model.createExpression(O_Sum);


                    // sum over programs
                    for (int P = 0; P < num_spots; P++) 
                    {
                        // building summand for P sum 
                        LSExpression term_logPU = model.createExpression(O_Prod, u[P], lsdouble(log_p_data(i,P)));
                        sum_logPU.addOperand(term_logPU);
                    }

                    // perform exponentation
                    LSExpression exp_logPU = model.createExpression(O_Exp, sum_logPU);
                    
                    // multiply by w_i
                    LSExpression term_wlogPU = model.createExpression(O_Prod, exp_logPU, lsdouble(w_data(i,0)));
                    sum_Individuals.addOperand(term_wlogPU);
                }

                // devide by target universe
                sum_Individuals = model.createExpression(O_Div, sum_Individuals, lsdouble(target_universe));

                // subtract summation from 1.0
                LSExpression Objective_function = model.createExpression(O_Sub, lsdouble(1.0), sum_Individuals);

                model.addObjective(Objective_function, OD_Maximize);
                model.close();

                // Local Solver settings
                LSPhase phase = localsolver.createPhase();
                phase.setTimeLimit(time_limit);

                // set the annealing level
                localsolver.getParam().setAnnealingLevel(annealing_level);
                
                std::cout << "\ncalling local solver" << std::endl;
                localsolver.solve();

                //solution.clear();
                std::cout << "\nprinting and returning solution" << std::endl;
                for (int P = 0; P < num_spots; ++P)
                {
                    solution[P] = static_cast<int>(u[P].getValue());
                }

            } catch (LSException *e) {
                cout << "LSException:" << e->getMessage() << std::endl;
                exit(1);
            }
        }

        void printSolution() 
        {
            for (unsigned int P = 0; P < num_spots; ++P)
            {
                std::cout << bp::extract<int>(solution[P]) << " ";
            }
            std::cout << std::endl;
        }

};


// Expose to Python
np::ndarray optimizeTRP
(
    np::ndarray const & grps,
    np::ndarray const & trps, 
    double const grp_limit, 
    int const time_limit = 10,
    int const annealing_level = 0
) 
{
    std::cout << "\nSet up local solver to Maximize TRPs" << std::endl;
    OptimizeTRP model(grps, trps, grp_limit);
    std::cout << "\nsolving" << std::endl;
    model.solve(time_limit, annealing_level);
    return model.solution;
}

np::ndarray optimizeReach
(
    np::ndarray const & log_p,
    np::ndarray const & w,
    np::ndarray const & grps, 
    double const grp_limit, 
    int const time_limit = 10,
    int const annealing_level = 0
) 
{
    std::cout << "\nSet up local solver to maximize Reach" << std::endl;
    OptimizeReach model(log_p, w, grps, grp_limit);
    std::cout << "\nsolving" << std::endl;
    model.solve(time_limit, annealing_level);
    return model.solution;
}


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
    bp::def("optimizeTRP", optimizeTRP);
    bp::def("optimizeReach", optimizeReach);
    bp::def("test", test_numpy_wrapper);
}
