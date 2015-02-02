#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include "localsolver.h"
#include "NumPyArrayData.hpp" 
#include <boost/numpy.hpp>
#include <boost/scoped_array.hpp>
#include <cmath>

using namespace localsolver;
using namespace std;
namespace bp = boost::python;
namespace np = boost::numpy;

class OptimizeGeneral {
    public:

        // constructor
        OptimizeGeneral
        (
            boost::python::dict const & prob_data_pydict,           // python dictionary of viewership probability
            boost::python::dict const & weight_data_pydict,         // python dictionary of weight data
            double              const   reach_power_py,             // strengh ot he reach term of the objective function 
            boost::python::dict const & constraints_data_pydict,    // python dictionary of constraints data
            boost::python::dict const & priorities_data_pydict,     // python dictionary of priority data, keyed by media_plan_id
            np::ndarray         const & np_u0,                      // Initial Conditions
            np::ndarray         const & np_recency_score,           // recency score
            double              const   np_recency_power,           // recency power
            np::ndarray         const & np_creative_duration_array, // integer array (0 or 1) of spot length data 
            np::ndarray         const & np_durations_histogram,     // desired normailized histograms of creative durations
            double              const   np_duration_power,          // recency power 
            np::ndarray         const & np_core_demo_TRPs           // recency score
        )
            : prob_data( prob_data_pydict)                // extract the python dict to a c++ boost::python::dict
            , weight_data( weight_data_pydict)            // extract the python dict to a c++ boost::python::dict
            , reach_power( reach_power_py)                           
            , constraints_data( constraints_data_pydict)  // extract the python dict to a c++ boost::python::dict
            , priorities_data( priorities_data_pydict)    // extract the python dict to a c++ boost::python::dict
            , media_plans( prob_data.keys())                              
            , u0( np_u0)
            , recency_score(np_recency_score)
            , recency_power(np_recency_power)
            , creative_duration_array(np_creative_duration_array)
            , durations_histogram( np_durations_histogram)
            , duration_power(np_duration_power)
            , core_demo_TRPs             (np_core_demo_TRPs) 
            , constraint_lhs_index       (bp::extract<np::ndarray>(constraints_data["lhs_index"]))
            , constraint_lhs_variable_id (bp::extract<np::ndarray>(constraints_data["lhs_variable_id"]))      
            , constraint_lhs_coeff       (bp::extract<np::ndarray>(constraints_data["lhs_coeff"]))
            , constraint_rhs             (bp::extract<np::ndarray>(constraints_data["rhs"]))
            , constraint_sense           (bp::extract<np::ndarray>(constraints_data["sense"]))
        {
            std::cout << "instantating Local Solver: Optimize General"  << std::endl;
            num_media_plans = bp::len(media_plans);

            // buld number of variables
            num_variables = 0;
            for(int i = 0; i < num_media_plans; i++)
            {
                // extract the probability data array for the given media plan
                np::ndarray const p_array = bp::extract<np::ndarray>(prob_data[media_plans[i]]);
                
                // get the number of columns (i.e. number of decision varibles in this media plan)
                int num_variables_mp = p_array.shape(1); 
                num_variables += num_variables_mp; 
            }
            std::cout << "number of variables: " << num_variables << std::endl;
        }

        /* Number of decision variables. */
        int num_variables;
        int num_media_plans;

        // data for reach  
        boost::python::dict const prob_data;
        boost::python::dict const weight_data;
        double const reach_power; 

        // constraints
        boost::python::dict const constraints_data;
        
        // priorities of media plan
        boost::python::dict const priorities_data;
        boost::python::list const media_plans;    // keys of the prob_data and weight_data dictionary

        // initial condtions
        np::ndarray u0;

        // initial condtions
        np::ndarray recency_score;
        double recency_power; 
        
        // creative duration data
        np::ndarray const creative_duration_array;  // 2D array, row index for media plan, column for duration
        np::ndarray const durations_histogram;     // 2D array, row index for media plan, column for duration
        double duration_power; 
        np::ndarray const core_demo_TRPs;           // core demo TRP
        
        // constraint data
        // LHS
        np::ndarray constraint_lhs_index; 
        np::ndarray constraint_lhs_variable_id; 
        np::ndarray constraint_lhs_coeff;      

        // RHS
        np::ndarray constraint_rhs; 

        // sense
        np::ndarray constraint_sense;

        /* LocalSolver. */
        LocalSolver localsolver;

        /* Decision variables, i.e. units awarded */
        std::vector<localsolver::LSExpression> u;

        // create decision variables u[i]
        void create_decision_variables()
        {
            LSModel model = localsolver.getModel();
            
            std::cout << "number of variables: " << num_variables << std::endl;
            u.resize(num_variables);

            // get variable type and extract upper and lower bounds
            for (int i = 0; i < num_variables; i++) 
            {
                u[i] = model.createExpression(O_Bool);

                stringstream s;
                s << "u[" << i << "]";
                u[i].setName(s.str());
            }
        }


        // Build Constraints
        void build_constraints()
        {

            LSModel model = localsolver.getModel();
            
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
                LSExpression current_constraint_LHS = model.createExpression(O_Sum);

                // extract LHS data and build LHS expression
                while( bp::extract<int>(constraint_lhs_index[sparce_matrix_row]) == c )
                {
                    int j        = bp::extract<int>(   constraint_lhs_variable_id[sparce_matrix_row]);
                    double value = bp::extract<double>(constraint_lhs_coeff[sparce_matrix_row]);
                    
                    // build and add the next term
                    LSExpression current_term_LHS = model.createExpression(O_Prod, u[j], lsdouble(value));
                    current_constraint_LHS.addOperand(current_term_LHS);

                    // increment the sparce matrix row
                    sparce_matrix_row++;

                    // if we are at the last row, break out of the while loop
                    if( sparce_matrix_row == num_sparce_matrix_rows)
                    {
                        break;
                    }
                }

                // extract the RHS 
                lsdouble rhs_value = lsdouble(bp::extract<double>(constraint_rhs[c]));

                // extract the sense
                int sense_int = bp::extract<int>(constraint_sense[c]);
                LSOperator sense = LSOperator(O_Eq);

                switch (sense_int)
                {
                  case -1:
                     sense = LSOperator(O_Leq);
                     break;
                  case 0:
                     sense = LSOperator(O_Eq);
                     break;
                  case 1:
                     sense = LSOperator(O_Geq);
                     break;
                  default:
                     break;
                }

                // add constraint to the model
                model.addConstraint(model.createExpression(sense, current_constraint_LHS, rhs_value));
            }
        }
       

       // build TRP objective function
       LSExpression build_objective_function_TRP()
       {
            LSModel model = localsolver.getModel();

            LSExpression TRP_objective = model.createExpression(O_Sum);

            int u_idx = 0;  // for indexing into the u[] vector
            // loop over media plans
            for (int mp = 0; mp < num_media_plans; mp++)
            {

                // extract the probability data array for the given media plan
                np::ndarray const p_array= bp::extract<np::ndarray>(prob_data[media_plans[mp]]);
                np::ndarray const w_array= bp::extract<np::ndarray>(weight_data[media_plans[mp]]);

                // use wrapper for easy indexing
                NumPyArrayData<double> p(p_array);
                NumPyArrayData<double> w(w_array);
            
                // get the number of columns (i.e. number of decision varibles in this media plan)
                int num_individuals_in_media_plan = p_array.shape(0); 
                int num_variables_in_media_plan   = p_array.shape(1); 

                // loop over decision variables in the current media plan
                for (int P = 0; P < num_variables_in_media_plan; P++)
                {
                    // compute the target average audience for variable P 
                    double targetAA_P  =  0.0; 
                    for (int i = 0; i < num_individuals_in_media_plan; i++)
                    {
                        targetAA_P += w(i,0)*p(i,P);
                    }
                    
                    // build LSEpression and add to the expression
                    LSExpression current_term = model.createExpression(O_Prod, u[u_idx], lsdouble(targetAA_P));
                    TRP_objective.addOperand(current_term);
                    u_idx++;
                }
            }

            return(TRP_objective); 
        }


        // build linear expansion objective function
        LSExpression build_objective_function_linear()
        {
            LSModel model = localsolver.getModel();

            LSExpression linear_objective = model.createExpression(O_Sum);

            int u_idx = 0;  // for indexing into the u[] vector
            // loop over media plans
            for (int mp = 0; mp < num_media_plans; mp++)
            {
                // extract the probability data array for the given media plan
                np::ndarray const p_array= bp::extract<np::ndarray>(prob_data[media_plans[mp]]);
                np::ndarray const w_array= bp::extract<np::ndarray>(weight_data[media_plans[mp]]);

                // use wrapper for easy indexing
                NumPyArrayData<double> p(p_array);
                NumPyArrayData<double> w(w_array);
            
                // get the number of columns (i.e. number of decision varibles in this media plan)
                int num_individuals_in_media_plan = p_array.shape(0); 
                int num_variables_in_media_plan   = p_array.shape(1); 

                // loop over decision variables in the current media plan
                for (int P = 0; P < num_variables_in_media_plan; P++)
                {
                    // compute the target average audience for variable P 
                    double targetAA_P  =  0.0; 
                    for (int i = 0; i < num_individuals_in_media_plan; i++)
                    {
                        double p_iP = p(i,P);
                        // make sure it's not 1.0
                        if(p_iP >= 1.0){
                            p_iP = .999;
                        }
                        targetAA_P += -1.0*w(i,0)*log(1.0 - p_iP);
                    }
                    
                    // build LSEpression and add to the expression
                    LSExpression current_term = model.createExpression(O_Prod, u[u_idx], lsdouble(targetAA_P));
                    linear_objective.addOperand(current_term);
                    u_idx++;
                }
            }

            return( linear_objective);
        }
        

        // build reach objective component
        LSExpression build_reach_expression()
        {
            std::cout << "building reach component of objective function"  << std::endl;
            
            int num_media_plans = bp::len(media_plans);

            LSModel model = localsolver.getModel();

            LSExpression reach_objective = model.createExpression(O_Sum);

            // keep track of shift in index for u[]
            int u_index_shift = 0;  
            
            // loop over media plans
            for (int mp = 0; mp < num_media_plans; mp++)
            {
                // media plan id
                int media_plan_id = bp::extract<int>(media_plans[mp]);

                // extract the probability data array for the given media plan
                np::ndarray const p_array= bp::extract<np::ndarray>(prob_data[media_plans[mp]]);
                np::ndarray const w_array= bp::extract<np::ndarray>(weight_data[media_plans[mp]]);

                // use wrapper for easy indexing
                NumPyArrayData<double> p(p_array);
                NumPyArrayData<double> w(w_array);
            
                // get the number of columns (i.e. number of decision varibles in this media plan)
                int num_individuals_in_media_plan = p_array.shape(0); 
                int num_variables_in_media_plan   = p_array.shape(1); 
               
                // build the possible universe for the given media plan
                double possible_universe = 0.0; 

                // build outer summation over individuals
                LSExpression sum_Individuals = model.createExpression(O_Sum);
                for (int i = 0; i < num_individuals_in_media_plan; i++)
                {

                    // update universe
                    possible_universe += w(i,0);

                    // for each individual, loop over the programs P
                    // build argument of the exponential
                    LSExpression sum_logPU = model.createExpression(O_Sum);
                    for (int P = 0; P < num_variables_in_media_plan; P++)
                    {
                        double p_iP = p(i,P);
                        // make sure it's not 1.0
                        if(p_iP >= 1.0){
                            p_iP = .999;
                        }

                        // building summand for P sum 
                        LSExpression term_logPU = model.createExpression(O_Prod, u[P + u_index_shift], lsdouble(log(1.0 - p_iP)));
                        sum_logPU.addOperand(term_logPU);
                    }
                    
                    // perform exponentation
                    LSExpression exp_logPU = model.createExpression(O_Exp, sum_logPU);
                    
                    // multiply by w_i
                    LSExpression term_wlogPU = model.createExpression(O_Prod, exp_logPU, lsdouble(w(i,0)));
                    sum_Individuals.addOperand(term_wlogPU);

                }

                
                // divide by target universe
                sum_Individuals = model.createExpression(O_Div, sum_Individuals, lsdouble(possible_universe));

                // subtract summation from 1.0
                LSExpression reach_per_media_plan = model.createExpression(O_Sub, lsdouble(1.0), sum_Individuals);

                // weight the reach by priority
                double priority = bp::extract<double>(priorities_data[media_plans[mp]]);
                reach_per_media_plan = model.createExpression(O_Prod, reach_per_media_plan, lsdouble(priority));
                
                // print for checking
                std::cout << "media_plan_id : " << media_plan_id << std::endl;
                std::cout << "\t" << possible_universe <<  " possible_universe" << std::endl;
                std::cout << "\tpriority = " << priority << std::endl;

                // add to total reach
                reach_objective.addOperand(reach_per_media_plan); 
                
                // update u_index_shift for next block-matrix
                u_index_shift += num_variables_in_media_plan;
            }
            
            std::cout << "done\n\n" << std::endl;

            reach_objective = model.createExpression(O_Div, reach_objective, lsdouble(num_media_plans)); 
            return( reach_objective);
        }

        // build recency objective component
        LSExpression build_recency_expression()
        {
            std::cout << "building recency component of objective function"  << std::endl;
            std::cout << "\trecency power = " << recency_power <<  std::endl;
            
            int num_media_plans = bp::len(media_plans);

            LSModel model = localsolver.getModel();

            LSExpression recency_objective = model.createExpression(O_Sum);


            // keep track of shift in index for u[]
            int u_index_shift = 0;  
            
            // loop over media plans
            for (int mp = 0; mp < num_media_plans; mp++)
            {
                // media plan id
                int media_plan_id = bp::extract<int>(media_plans[mp]);

                // extract the probability data array for the given media plan
                np::ndarray const p_array= bp::extract<np::ndarray>(prob_data[media_plans[mp]]);

                // extract the probability data array for the given media plan
                int num_variables_in_media_plan   = p_array.shape(1); 

                // for each media plan, compute reency score
                LSExpression recency_per_media_plan = model.createExpression(O_Sum);

                // set up for loop
                int P_start = u_index_shift;
                int P_stop  = P_start + num_variables_in_media_plan;

                std::cout << "\tP_start: " << P_start;
                std::cout << "\tP_stop: "  << P_stop  << std::endl;

                // loop over decision variables
                for (int P = P_start; P < P_stop; P++)
                {

                    // extract the recency score, times priority, for each decision variable
                    double r_P = bp::extract<double>(recency_score[P]);
                    LSExpression current_term = model.createExpression(O_Prod, u[P], lsdouble(r_P));
                    recency_per_media_plan.addOperand(current_term);
                }

                // mulitply by priority!
                double priority = bp::extract<double>(priorities_data[media_plans[mp]]);
                recency_per_media_plan = model.createExpression(O_Prod, recency_per_media_plan, lsdouble(priority)); 

                // print for checking
                std::cout << "media_plan_id : " << media_plan_id << std::endl;
                std::cout << "\tpriority = " << priority << std::endl;
                
                // add the reach to the total reach, after weighting by priority
                recency_objective.addOperand(recency_per_media_plan); 
                
                // update u_index_shift for next block-matrix
                u_index_shift += num_variables_in_media_plan;
            }

            std::cout << "done\n\n" << std::endl;

            recency_objective = model.createExpression(O_Div, recency_objective, lsdouble(num_media_plans)); 
            return(recency_objective);
        }

        // build creative duration expression
        LSExpression build_creative_duration_expression()
        {
            std::cout << "building creative_duration component of objective function"  << std::endl;
            std::cout << "\tduration power = " << duration_power <<  std::endl;

            LSModel model = localsolver.getModel();

            LSExpression creative_duration_objective = model.createExpression(O_Sum);
            
            // use wrapper for easy indexing
            NumPyArrayData<long> f_array(creative_duration_array);      // f data tags the creative duration
            NumPyArrayData<double> n_desired(durations_histogram);    // n_true is the desired histogram

            // helper variaables
            size_t num_media_plans = bp::len(media_plans);
            size_t num_durations   = creative_duration_array.shape(1);

            // keep track of shift in index for u[]
            int u_index_shift = 0;  
            
            // loop over media plans
            for (int mp = 0; mp < num_media_plans; mp++)
            {
                // get the media plan's priority 
                double priority = bp::extract<double>(priorities_data[media_plans[mp]]);

                // extract the probability data array for the given media plan
                np::ndarray const p_array= bp::extract<np::ndarray>(prob_data[media_plans[mp]]);

                // extract the probability data array for the given media plan
                int num_variables_in_media_plan   = p_array.shape(1); 

                int P_start = u_index_shift;
                int P_stop  = P_start + num_variables_in_media_plan;

                std::cout << "\tP_start: " << P_start;
                std::cout << "\tP_stop: "  << P_stop  << std::endl;

                // store the "numerators"
                std::vector<LSExpression> num;
                // initialize the std::vector 
                for (int d = 0; d < num_durations; d++)
                {
                    num.push_back(model.createExpression(O_Sum));
                }

                // store the "denominator
                LSExpression den = model.createExpression(O_Sum);
                den.addOperand(model.createConstant(lsdouble(1.0e-9)));
                
                // loop over decision variables
                for (int P = P_start; P < P_stop; P++)
                {
                    // extract the core_demo_AA
                    double core_demo_TRP = bp::extract<double>(core_demo_TRPs[P]);
                    // update denominator
                    den.addOperand(model.createExpression(O_Prod,u[P],lsdouble(core_demo_TRP)));


                    // for each duration, update numerator
                    for (int d = 0; d < num_durations; d++)
                    {
                        num[d].addOperand(model.createExpression(O_Prod,u[P],lsdouble(core_demo_TRP*f_array(P,d))));
                    }
                }
                // for each media plan, compute reency score
                 LSExpression score_per_media_plan = model.createExpression(O_Sum);
                
                // now that the histogram is built, construct store
                for (int d = 0; d < num_durations; d++)
                {
                    LSExpression n_hat       = model.createExpression(O_Div,num[d], den);
                    LSExpression n_diff      = model.createExpression(O_Sub, n_hat, lsdouble(n_desired(mp, d)));
                    LSExpression abs_n_diff  = model.createExpression(O_Abs, n_diff); 
                    abs_n_diff  = model.createExpression(O_Sqrt, abs_n_diff); 

                    //LSExpression abs_n_diff  = model.createExpression(O_Pow, n_diff, lsint(2)); 

                    //score_per_media_plan.addOperand(abs_n_diff);
                    score_per_media_plan.addOperand(abs_n_diff);
                }

                // make positive
                score_per_media_plan = model.createExpression(O_Prod,score_per_media_plan, lsdouble(-0.25));
                score_per_media_plan = model.createExpression(O_Sum, model.createConstant(lsdouble(1.0)), score_per_media_plan);

                // add the score for each media plan to the full traffic plan solution
                // scaled by priority
                score_per_media_plan = model.createExpression(O_Prod,score_per_media_plan, lsdouble(priority));
                score_per_media_plan = model.createExpression(O_Div,score_per_media_plan, lsdouble(num_media_plans));

                // add the media plan score the totals
                creative_duration_objective.addOperand(score_per_media_plan); 
                
                // update u_index_shift for next block-matrix
                u_index_shift += num_variables_in_media_plan;
            }

            std::cout << "done\n\n" << std::endl;
            return(creative_duration_objective);
        }
        

        // solver the model
        void solve(int time_limit, int annealing_level) {
            try 
            {
                LSModel model = localsolver.getModel();

                create_decision_variables();

                build_constraints();
               
                //LSExpression TRP_objective =  build_objective_function_TRP();
                //LSExpression linear_objective =  build_objective_function_linear();

                // build objective function
                LSExpression reach_term = build_reach_expression();
                //LSExpression reach_term = build_equivalized_reach_expression();
                reach_term =  model.createExpression(O_Prod, reach_term, lsdouble(reach_power));

                // recency
                LSExpression recency_term = build_recency_expression();
                recency_term =  model.createExpression(O_Prod, recency_term, lsdouble(recency_power));

                // creative duration
                LSExpression creative_duration_term = build_creative_duration_expression();
                creative_duration_term =  model.createExpression(O_Prod,creative_duration_term, lsdouble(duration_power));

                // add the terms together
                LSExpression nonlinear_objective = model.createExpression(O_Sum, reach_term, 
                                                                                 recency_term,
                                                                                 creative_duration_term);

                // commit the objective function
                model.addObjective(nonlinear_objective, OD_Maximize);
                //model.addObjective(creative_duration_term, OD_Maximize);
                model.close();

                // Local Solver settings
                LSPhase phase = localsolver.createPhase();
                phase.setTimeLimit(time_limit);

                // set the annealing level
                localsolver.getParam().setAnnealingLevel(annealing_level);

                // set initial conditions
                for (int P = 0; P < num_variables; P++) 
                {
                    NumPyArrayData<long> u0_data(u0);
                    u[P].setValue(lsint(u0_data(P)));
                }

                
                std::cout << "\ncalling local solver" << std::endl;
                localsolver.solve();
                
                //std::cout << "creative_duration_term = " << creative_duration_term.getValue() << std::endl;

            } catch (LSException *e) {
                cout << "LSException:" << e->getMessage() << std::endl;
                exit(1);
            }
        }

        np::ndarray getSolution() 
        {
        
            np::ndarray solution = np::zeros(bp::make_tuple(num_variables),  np::dtype::get_builtin<int>());
            
            for (int P = 0; P < num_variables; ++P)
            {
                solution[P] = static_cast<int>(u[P].getValue());
            }

            return solution;
        }

};

class OptimizeLinear {
    public:

        // constructor
        OptimizeLinear
        (
            np::ndarray const & np_obj_coeff,                    // the linear coefficients of the objective function
            boost::python::dict const & constraints_data_pydict, // python dictionary of constraints data
            np::ndarray const & np_u0                            // Initial Conditions
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

        /* LocalSolver. */
        LocalSolver localsolver;

        /* Decision variables, i.e. units awarded */
        std::vector<localsolver::LSExpression> u;

        // create decision variables u[i]
        void create_decision_variables()
        {
            LSModel model = localsolver.getModel();
            
            std::cout << "number of variables: " << num_variables << std::endl;
            u.resize(num_variables);

            // get variable type and extract upper and lower bounds
            for (int i = 0; i < num_variables; i++) 
            {
                u[i] = model.createExpression(O_Bool);

                stringstream s;
                s << "u[" << i << "]";
                u[i].setName(s.str());
            }
        }


        // Build Constraints
        void build_constraints()
        {

            LSModel model = localsolver.getModel();
            
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
                LSExpression current_constraint_LHS = model.createExpression(O_Sum);

                // extract LHS data and build LHS expression
                while( bp::extract<int>(constraint_lhs_index[sparce_matrix_row]) == c )
                {
                    int j        = bp::extract<int>(   constraint_lhs_variable_id[sparce_matrix_row]);
                    double value = bp::extract<double>(constraint_lhs_coeff[sparce_matrix_row]);
                    
                    // build and add the next term
                    LSExpression current_term_LHS = model.createExpression(O_Prod, u[j], lsdouble(value));
                    current_constraint_LHS.addOperand(current_term_LHS);

                    // increment the sparce matrix row
                    sparce_matrix_row++;

                    // if we are at the last row, break out of the while loop
                    if( sparce_matrix_row == num_sparce_matrix_rows)
                    {
                        break;
                    }
                }

                // extract the RHS 
                lsdouble rhs_value = lsdouble(bp::extract<double>(constraint_rhs[c]));

                // extract the sense
                int sense_int = bp::extract<int>(constraint_sense[c]);
                LSOperator sense = LSOperator(O_Eq);

                switch (sense_int)
                {
                  case -1:
                     sense = LSOperator(O_Leq);
                     break;
                  case 0:
                     sense = LSOperator(O_Eq);
                     break;
                  case 1:
                     sense = LSOperator(O_Geq);
                     break;
                  default:
                     break;
                }

                // add constraint to the model
                model.addConstraint(model.createExpression(sense, current_constraint_LHS, rhs_value));
            }
        }
       

       // build TRP objective function
       LSExpression build_objective_function()
       {
            LSModel model = localsolver.getModel();

            LSExpression linear_objective = model.createExpression(O_Sum);


            // loop over decision variables in the current media plan
            for (int P = 0; P < num_variables; P++)
            {
                // build LSEpression and add to the expression
                LSExpression current_term = model.createExpression(O_Prod, u[P], lsdouble(bp::extract<double>(obj_coeff[P])));
                linear_objective.addOperand(current_term);
            }

            return(linear_objective); 
        }


        // solver the model
        void solve(int time_limit, int annealing_level) {
            try 
            {
                LSModel model = localsolver.getModel();

                create_decision_variables();

                // build constraints
                build_constraints();

                // objective function
                LSExpression linear_objective =  build_objective_function();

                // commit the objective function
                model.addObjective(linear_objective, OD_Maximize);
                model.close();

                // Local Solver settings
                LSPhase phase = localsolver.createPhase();
                phase.setTimeLimit(time_limit);

                // set the annealing level
                localsolver.getParam().setAnnealingLevel(annealing_level);

                // set initial conditions
                for (int P = 0; P < num_variables; P++) 
                {
                    NumPyArrayData<long> u0_data(u0);
                    u[P].setValue(lsint(u0_data(P)));
                }
                
                std::cout << "\ncalling local solver" << std::endl;
                localsolver.solve();

            } catch (LSException *e) {
                cout << "LSException:" << e->getMessage() << std::endl;
                exit(1);
            }
        }

        np::ndarray getSolution() 
        {
        
            np::ndarray solution = np::zeros(bp::make_tuple(num_variables),  np::dtype::get_builtin<int>());
            
            for (int P = 0; P < num_variables; ++P)
            {
                solution[P] = static_cast<int>(u[P].getValue());
            }

            return solution;
        }

};

// Expose to Python
np::ndarray optimizeGeneral
(
    boost::python::dict   const & prob_data_pydict,           // python dictionary of viewership probability
    boost::python::dict   const & weight_data_pydict,         // python dictionary of weight data
    double                const   reach_power,                // strengh ot he reach term of the objective function 
    boost::python::dict   const & constraints_data_pydict,    // python dictionary of constraints data
    boost::python::dict   const & priorities_data_pydict,     // python dictionary of priorities data
    np::ndarray           const & u0,                         // initial condtions
    np::ndarray           const & recency,                    // (math) column vector of recency scores
    double                const & recency_power,              // strength of the recency term of the objective function
    np::ndarray           const & creative_duration_array,    // integer array (0 or 1) of spot length data 
    np::ndarray           const & durations_histogram,        // desired normailized histograms of creative durations
    double                const   duration_power,             // strength of the creative_duration term in the objective function
    np::ndarray           const & core_demo_TRPs,             // core demo average audience, TRP 
    int                   const   time_limit = 10,
    int                   const   annealing_level = 0
) 
{
    std::cout << "\nSet up local solver to Maximize full reach" << std::endl;
    OptimizeGeneral model(prob_data_pydict, 
                          weight_data_pydict, 
                          reach_power, 
                          constraints_data_pydict,
                          priorities_data_pydict, 
                          u0,
                          recency, 
                          recency_power, 
                          creative_duration_array,
                          durations_histogram, 
                          duration_power, 
                          core_demo_TRPs
                          );
    std::cout << "\nsolving" << std::endl;
    model.solve(time_limit, annealing_level);
    return model.getSolution();
}

// Expose to Python
np::ndarray optimizeLinear
(
    np::ndarray           const & obj_coeff,                  // the linear coefficients of the objective function
    boost::python::dict   const & constraints_data_pydict,    // python dictionary of constraints data
    np::ndarray const & u0,
    int const time_limit = 10,
    int const annealing_level = 0
) 
{
    std::cout << "\nSet up local solver to Maximize Gurobi Linear Problems" << std::endl;
    OptimizeLinear model(obj_coeff, constraints_data_pydict, u0);
    std::cout << "\nsolving" << std::endl;
    model.solve(time_limit, annealing_level);
    return model.getSolution();
}


BOOST_PYTHON_MODULE(localsolver_tools) {
    np::initialize();  // have to put this in any module that uses Boost.NumPy
    bp::def("optimizeLinear", optimizeLinear);
    bp::def("optimizeGeneral", optimizeGeneral);
}
