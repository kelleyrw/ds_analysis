#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include "NumPyArrayData.hpp" 
#include <boost/numpy.hpp>
#include <boost/scoped_array.hpp>
#include <cmath>
#include "gurobi_c++.h"

using namespace std;
namespace bp = boost::python;
namespace np = boost::numpy;

class OptimizeGurobi {
    public:

        // constructor
        OptimizeGurobi
        (
            np::ndarray const         & np_obj_coeff,            // the linear coefficients of the objective function
            boost::python::dict const & constraints_data_pydict, // python dictionary of constraints data
            np::ndarray const         & np_u0                    // Initial Conditions
        )
            : num_variables(np_obj_coeff.shape(0))
            , obj_coeff(np_obj_coeff) 
            , u0(np_u0)
            , constraint_lhs_index       (bp::extract<np::ndarray>(constraints_data_pydict["lhs_index"]))
            , constraint_lhs_variable_id (bp::extract<np::ndarray>(constraints_data_pydict["lhs_variable_id"]))      
            , constraint_lhs_coeff       (bp::extract<np::ndarray>(constraints_data_pydict["lhs_coeff"]))
            , constraint_rhs             (bp::extract<np::ndarray>(constraints_data_pydict["rhs"]))
            , constraint_sense           (bp::extract<np::ndarray>(constraints_data_pydict["sense"]))
        {
            std::cout << "instantating Local Solver: Optimize Linear"  << std::endl;
            std::cout << "number of variables: " << num_variables << std::endl;
        }

        /* Number of decision variables. */
        int num_variables;
        int num_media_plans;

        // objective function data 
        np::ndarray obj_coeff;
        
        // boost::python::dict const constraints_data;
        boost::python::list const media_plans;    // keys of the prob_data and weight_data dictionary

        // initial condtions
        np::ndarray u0;

        // constraint data
        // LHS
        np::ndarray constraint_lhs_index; 
        np::ndarray constraint_lhs_variable_id; 
        np::ndarray constraint_lhs_coeff;      

        // RHS
        np::ndarray constraint_rhs; 

        // sense
        np::ndarray constraint_sense;

        /* gurobi. */
        GRBModel model; 
        
        /* Decision variables, i.e. units awarded */
        std::vector<GRBVar> u;

        // create decision variables u[i]
        void create_decision_variables()
        {
            
            std::cout << "number of variables: " << num_variables << std::endl;
            u.resize(num_variables);

            // get variable type and extract upper and lower bounds
            for (int P = 0; P < num_variables; P++) 
            {
                stringstream s;
                s << "u[" << P << "]";
                GRBVar u[P] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, s);
            }
        }


        // Build Constraints
        void build_constraints()
        {

            std::cout << "\nbuilding constraints\n" << std::endl;
            int num_sparce_matrix_rows = constraint_lhs_coeff.shape(0);
            int num_constraints = constraint_rhs.shape(0);

            // ndarray accessor wrappers
            // only use for 2D and higher dimensional arrays
            NumPyArrayData<double> constraint_lhs_coeff_data(constraint_lhs_coeff);

            // point to current row in sparce matrix 
            int sparce_matrix_row = 0;

            // loop over constraints
            for(int c = 0; c < num_constraints; c++)
            {
                //cout << "\tadding constraint " << c <<  endl;
                // initialize LHS expressions
                GRBLinExpr	current_constraint_LHS = 0.0;

                // extract LHS data and build LHS expression
                while( bp::extract<int>(constraint_lhs_index[sparce_matrix_row]) == c )
                {
                    int j        = bp::extract<int>(   constraint_lhs_variable_id[sparce_matrix_row]);
                    double value = bp::extract<double>(constraint_lhs_coeff[sparce_matrix_row]);
                    
                    // build and add the next term
                    current_constraint_LHS += u[j]*value;

                    // increment the sparce matrix row
                    sparce_matrix_row++;

                    // if we are at the last row, break out of the while loop
                    if( sparce_matrix_row == num_sparce_matrix_rows)
                    {
                        break;
                    }
                }

                // extract the RHS 
                GRBLinExpr rhs_value = bp::extract<double>(constraint_rhs[c]);

                // extract the sense
                int sense_int = bp::extract<int>(constraint_sense[c]);
                std::string sense;

                switch (sense_int)
                {
                  case -1:
                     sense = "GRB_LESS_EQUAL";
                     break;
                  case 0:
                     sense = "GRB_EQUAL";
                     break;
                  case 1:
                     sense = "GRB_GREATER_EQUAL";
                     break;
                  default:
                     break;
                }

                // add constraint to the model
                model.addConstr(	current_constraint_LHS, sense.c_str(), rhs_value);
            }
        }

       // build objective function
       GRBLinExpr build_objective_function()
       {
            GRBLinExpr linear_objective = 0.0;

            // loop over decision variables in the current media plan
            for (int P = 0; P < num_variables; P++)
            {
                linear_objective += lsdouble(bp::extract<double>(obj_coeff[P])) * u[P];
            }

            return(linear_objective); 
        }


        // solver the model
        void solve(int time_limit, int annealing_level) {
            try 
            {

                create_decision_variables();
                build_constraints();

                // objective function
                linear_objective =  build_objective_function();
                model.setObjective(linear_objective, GRB_MINIMIZE);

                // commit the objective function
                model.update;

                // // set initial conditions
                // for (int P = 0; P < num_variables; P++) 
                // {
                //     NumPyArrayData<long> u0_data(u0);
                //     u[P].setValue(lsint(u0_data(P)));
                // }
                
                std::cout << "\noptimizing in Gurobi" << std::endl;
                model.optimize();

            } catch(GRBException e) {
              cout << "Error code = " << e.getErrorCode() << endl;
              cout << e.getMessage() << endl;
            } catch(...) {
              cout << "Exception during optimization" << endl;
            }
        }

        np::ndarray getSolution() 
        {
            np::ndarray solution = np::zeros(bp::make_tuple(num_variables),  np::dtype::get_builtin<int>());
            
            for (int P = 0; P < num_variables; ++P)
            {
                solution[P] = static_cast<int>(u[P].get(GRB_DoubleAttr_X));
            }

            return solution;
        }

};

// Expose to Python
np::ndarray optimizeGurobi
(
    np::ndarray           const & obj_coeff,                  // the linear coefficients of the objective function
    boost::python::dict   const & constraints_data_pydict,    // python dictionary of constraints data
    np::ndarray const & u0,
    int const time_limit = 10
) 
{
    std::cout << "\nSet up Gurobi for optimiztion" << std::endl;
    OptimizeGurobi model(obj_coeff, constraints_data_pydict, u0);
    std::cout << "\nsolving" << std::endl;
    model.solve(time_limit);
    return model.getSolution();
}

BOOST_PYTHON_MODULE(gurobi_tools) {
    np::initialize();  // have to put this in any module that uses Boost.NumPy
    bp::def("optimizeGurobi", optimizeLinear);
}
