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
            np::ndarray          const & np_obj_coeff,            // the linear coefficients of the objective function
            np::ndarray          const & np_lb,                   // lower bound
            np::ndarray          const & np_ub,                   // upper bound
            np::ndarray          const & np_vtype,                // variable type
            boost::python::dict  const & constraints_data_pydict, // python dictionary of constraints data
            np::ndarray          const & np_u0,                   // Initial Conditions
            int                  const & bp_model_sense, 
            bool                 const & use_compute_server,       // Gurobi Compute Server
            double               const & time_limit
        )
            : model_sense (bp_model_sense)
            , obj_coeff(np_obj_coeff) 
            , lb(np_lb) 
            , ub(np_ub) 
            , vtype(np_vtype) 
            , constraint_lhs_index       (bp::extract<np::ndarray>(constraints_data_pydict["lhs_index"]))
            , constraint_lhs_variable_id (bp::extract<np::ndarray>(constraints_data_pydict["lhs_variable_id"]))      
            , constraint_lhs_coeff       (bp::extract<np::ndarray>(constraints_data_pydict["lhs_coeff"]))
            , constraint_rhs             (bp::extract<np::ndarray>(constraints_data_pydict["rhs"]))
            , constraint_sense           (bp::extract<np::ndarray>(constraints_data_pydict["sense"]))
            , num_variables(np_obj_coeff.shape(0))
            , u0(np_u0)
            , env(((use_compute_server) ? 
                    GRBEnv("", "ec2-54-160-64-57.compute-1.amazonaws.com", -1, "", 0, 10) :
                    GRBEnv() ))
            , model(GRBModel(env))
        {
            std::cout << "instantating: Optimize Gurobi"  << std::endl;
            std::cout << "number of variables: " << num_variables << std::endl;
            
            // set the time limit
            model.getEnv().set(GRB_DoubleParam_TimeLimit, time_limit); 
        }

        // constructor
        OptimizeGurobi
        (
            boost::python::dict  const & objectives_data_pydict, // python dictionary of constraints data
            boost::python::dict  const & constraints_data_pydict, // python dictionary of constraints data
            np::ndarray          const & np_u0,                   // Initial Conditions
            bool                 const & use_compute_server,      // Gurobi Compute Server
            double               const & time_limit
        )
            : model_sense                (bp::extract<int>(objectives_data_pydict["model_sense"]))
            , obj_coeff                  (bp::extract<np::ndarray>(objectives_data_pydict["obj_coeff"]))
            , lb                         (bp::extract<np::ndarray>(objectives_data_pydict["lb"]))
            , ub                         (bp::extract<np::ndarray>(objectives_data_pydict["ub"]))
            , vtype                      (bp::extract<np::ndarray>(objectives_data_pydict["vtype"]))
            , constraint_lhs_index       (bp::extract<np::ndarray>(constraints_data_pydict["lhs_index"]))
            , constraint_lhs_variable_id (bp::extract<np::ndarray>(constraints_data_pydict["lhs_variable_id"]))      
            , constraint_lhs_coeff       (bp::extract<np::ndarray>(constraints_data_pydict["lhs_coeff"]))
            , constraint_rhs             (bp::extract<np::ndarray>(constraints_data_pydict["rhs"]))
            , constraint_sense           (bp::extract<np::ndarray>(constraints_data_pydict["sense"]))
            , num_variables              (obj_coeff.shape(0))
            , u0                         (np_u0)
            , env(((use_compute_server) ? 
                    GRBEnv("", "ec2-54-160-64-57.compute-1.amazonaws.com", -1, "", 0, 10) :
                    GRBEnv() ))
            , model(GRBModel(env))
        {
            std::cout << "instantating: Optimize Gurobi"  << std::endl;
            std::cout << "number of variables: " << num_variables << std::endl;
            
            // set the time limit
            model.getEnv().set(GRB_DoubleParam_TimeLimit, time_limit); 
        }

        // objective function data 
        const int model_sense;
        const np::ndarray obj_coeff;
        const np::ndarray lb;
        const np::ndarray ub;
        const np::ndarray vtype;

        // constraint data
        // LHS
        const np::ndarray constraint_lhs_index; 
        const np::ndarray constraint_lhs_variable_id; 
        const np::ndarray constraint_lhs_coeff;      
        const np::ndarray constraint_rhs; 
        const np::ndarray constraint_sense;

        /* Number of decision variables. */
        int num_variables;
        
        // initial condtions
        np::ndarray u0;

        /* gurobi. */
        GRBEnv env; 
        GRBModel model; //  = GRBModel(env); 
        
        /* Decision variables, i.e. units awarded */
        std::vector<GRBVar> u;

        // create decision variables 
        void create_decision_variables()
        {

            std::cout << "\ncreating decision variables" << std::endl;
            std::cout << "\ncreating objective function" << std::endl;
            u.resize(num_variables);

            // get variable type and extract upper and lower bounds
            for (int P = 0; P < num_variables; P++) 
            {
                // extract bounds
                double obj = bp::extract<double>(obj_coeff[P]);

                // extract bounds
                double lb_P = bp::extract<double>(lb[P]);
                double ub_P = bp::extract<double>(ub[P]);

                // extract the sense
                int vtype_int = bp::extract<int>(vtype[P]);
                std::string vtype_string;

                // string name
                stringstream s;
                s << "u[" << P << "]";

                switch (vtype_int)
                {
                  // BINARY
                  case 0: 
                     vtype_string = "GRB_BINARY";
                     u[P] = model.addVar(lb_P, ub_P, obj, GRB_BINARY, s.str());
                     break;

                  // INTEGER
                  case 1:
                     vtype_string = "GRB_INTEGER";
                     u[P] = model.addVar(lb_P, ub_P, obj, GRB_INTEGER, s.str());
                     break;

                  // CONTINUOUS
                  case 2:
                     vtype_string = "GRB_CONTINUOUS";
                     u[P] = model.addVar(lb_P, ub_P, obj, GRB_CONTINUOUS, s.str());
                     break;
                  default:
                     break;
                }

                // cout << "adding variable " << P <<  endl;
                // cout << "\tObj Coeff = "   << obj   << endl;
                // cout << "\tlb = "          << lb_P    << endl;
                // cout << "\tub = "          << ub_P    << endl;
                // cout << "\tvtype = "       << vtype_string << endl;
            }
            model.update();
        }


        // Build Constraints
        void build_constraints()
        {
            std::cout << "\nbuilding constraints" << std::endl;
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
                GRBLinExpr rhs_value = GRBLinExpr(bp::extract<double>(constraint_rhs[c]));

                // extract the sense
                int sense_int = bp::extract<int>(constraint_sense[c]);
                std::string sense;

                switch (sense_int)
                {
                  case -2:
                     sense = "GRB_LESS_EQUAL";
                     model.addConstr(current_constraint_LHS, GRB_LESS_EQUAL, rhs_value);
                     break;
                  case 0:
                     sense = "GRB_EQUAL";
                     model.addConstr(current_constraint_LHS, GRB_EQUAL, rhs_value);
                     break;
                  case 2:
                     sense = "GRB_GREATER_EQUAL";
                     model.addConstr(current_constraint_LHS, GRB_GREATER_EQUAL, rhs_value);
                     break;
                  default:
                     break;
                }
            }

            // update model
            model.update();
        }

        // solver the model
        void solve() {
            try 
            {
                std::cout << "\noptimizing....in C++" << std::endl;

                create_decision_variables();
                build_constraints();

                // sense of the optimization to maximize
                model.set(GRB_IntAttr_ModelSense, model_sense);
                model.update();

                // // set initial conditions
                // for (int P = 0; P < num_variables; P++) 
                // {
                //     NumPyArrayData<long> u0_data(u0);
                //     u[P].setValue(lsint(u0_data(P)));
                // }
                
                std::cout << "\nsending model to Gurobi compute server\n" << std::endl;
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
            // 
            // for (int P = 0; P < num_variables; ++P)
            // {
            //     solution[P] = static_cast<int>(u[P].get(GRB_DoubleAttr_X));
            // }
   
            // extract data
            double* array = model.get(GRB_DoubleAttr_X, model.getVars(), num_variables);
            std::vector<double> solution(array, array+num_variables);

            // convert to ndarray
            np::ndarray np_solution = np::zeros(bp::make_tuple(num_variables),  np::dtype::get_builtin<int>());
            std::copy(solution.begin(), solution.end(), reinterpret_cast<int*>(np_solution.get_data()));

            // delete point
            delete[] array;

            return np_solution;

        }

};

// Expose to Python
np::ndarray optimizeGurobi
(
    boost::python::dict   const & objective_data_pydict,   // python dictionary of constraints data
    boost::python::dict   const & constraints_data_pydict, // python dictionary of constraints data
    np::ndarray           const & u0,
    bool                  const   use_compute_server,    
    double                const   time_limit = 10.0
) 
{
    std::cout << "\nSet up Gurobi for optimiztion" << std::endl;
    OptimizeGurobi model(objective_data_pydict, 
                     constraints_data_pydict, 
                     u0, 
                     use_compute_server, 
                     time_limit
                     );
    model.solve();
    std::cout << "\ndone" << std::endl;
    return model.getSolution();
}

// Expose to Python
np::ndarray optimizeGurobi2
(
    np::ndarray           const & obj_coeff, 
    np::ndarray           const & lb,   
    np::ndarray           const & ub,    
    np::ndarray           const & vtype, 
    boost::python::dict   const & constraints_data_pydict, 
    int                   const   model_sense, 
    np::ndarray           const & u0,
    bool                  const   use_compute_server,    
    double                const   time_limit = 10.0
) 
{
    std::cout << "\nSet up Gurobi for optimiztion" << std::endl;
    OptimizeGurobi model(obj_coeff, 
                     lb,
                     ub,
                     vtype,
                     constraints_data_pydict, 
                     u0, 
                     model_sense, 
                     use_compute_server, 
                     time_limit
                     );
    model.solve();
    std::cout << "\ndone" << std::endl;
    return model.getSolution();
}

BOOST_PYTHON_MODULE(gurobi_tools) {
    np::initialize();  // have to put this in any module that uses Boost.NumPy
    bp::def("optimizeGurobi",  optimizeGurobi);
    bp::def("optimizeGurobi2", optimizeGurobi2);
}
