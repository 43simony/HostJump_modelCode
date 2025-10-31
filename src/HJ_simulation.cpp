//
//  HJ_simulation.cpp
//  
//
//  Created by Brandon Simony on 11/9/23.
//

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <chrono>
#include <random>
#include <sstream>
#include <cassert>
#include <algorithm>
#include "gsl/gsl_randist.h"

/*///////////////////////
//// GSL probability ////
////  Distributions  ////
///////////////////////*/

//For initiating and seeding gsl random number generator.
class Gsl_rng_wrapper
{
    gsl_rng* r;
    public:
        Gsl_rng_wrapper()
        {
            std::random_device rng_dev; //To generate a safe seed.
            long seed = time(NULL)*rng_dev();
            const gsl_rng_type* rng_type = gsl_rng_default;
            r = gsl_rng_alloc(rng_type);
            gsl_rng_set(r, seed);
        }
        ~Gsl_rng_wrapper() { gsl_rng_free(r); }
        gsl_rng* get_r() { return r; }
};


//Uniform RV using gsl.
double draw_uniform_gsl(double lo, double hi)
{
    static Gsl_rng_wrapper rng_w;
    static gsl_rng* r;
    if(r == nullptr){ r = rng_w.get_r(); }
    return gsl_ran_flat(r, lo, hi);
}


//Binomial RV using gsl.
int draw_binom_gsl(int N, double prob)
{
    static Gsl_rng_wrapper rng_w;
    static gsl_rng* r;
    if(r == nullptr){ r = rng_w.get_r(); }
    return gsl_ran_binomial(r, prob, N);
}


//NegativeBinomial RV using gsl.
int draw_negBinom_gsl(int N, double prob)
{
    static Gsl_rng_wrapper rng_w;
    static gsl_rng* r;
    if(r == nullptr){ r = rng_w.get_r(); }
    return gsl_ran_binomial(r, prob, N);
}


// Poisson RV using gsl. parameterized by mean
int draw_poisson_gsl(double lambda)
{
    static Gsl_rng_wrapper rng_w;
    static gsl_rng* r;
    if(r == nullptr){ r = rng_w.get_r(); }
    return gsl_ran_poisson(r, lambda);
}


//Beta RV using gsl.
double draw_beta_gsl(double a, double b)
{
    static Gsl_rng_wrapper rng_w;
    static gsl_rng* r;
    if(r == nullptr){ r = rng_w.get_r(); }
    return gsl_ran_beta(r, a, b);
}

//Beta RV using gsl.
double draw_beta_mix(std::vector<double> a, std::vector<double> b, std::vector<double> weights)
{
    static Gsl_rng_wrapper rng_w;
    static gsl_rng* r;
    if(r == nullptr){ r = rng_w.get_r(); }
    
    int idx = 0;
    double val = 0;
    double U1 = draw_uniform_gsl(0.0, 1.0);

    while(val == 0){
        if(U1 < weights[idx]){ val = draw_beta_gsl(a[idx], b[idx]); }
        idx++;
    }

    return val;
}


//Exponential RV using gsl. parameterized by a mean
double draw_exp_gsl(double mu)
{
    static Gsl_rng_wrapper rng_w;
    static gsl_rng* r;
    if(r == nullptr){ r = rng_w.get_r(); }
    return gsl_ran_exponential(r, mu);
}


//Gamma RV using gsl. k, theta parameterization
double draw_gamma_gsl(double shape, double scale)
{
    static Gsl_rng_wrapper rng_w;
    static gsl_rng* r;
    if(r == nullptr){ r = rng_w.get_r(); }
    return gsl_ran_gamma(r, shape, scale);
}


//vector-defined discrete probability distribution
double draw_discrete_gsl(std::vector<double> nodes, std::vector<double> weights)
{
    // ensure lengths of nodes and weights match
    if(nodes.size() != weights.size()){
        std::cout << "warning: node and weight vectors have different lengths" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    static Gsl_rng_wrapper rng_w;
    static gsl_rng* r;
    if(r == nullptr){ r = rng_w.get_r(); }
    double P1 = gsl_ran_flat(r, 0.0, 1.0);
    
    // Find the index where the random value falls in the cumulative distribution
    auto it = std::lower_bound(weights.begin(), weights.end(), P1); // returns weight value
    int index = std::distance(weights.begin(), it); // return position of the selected weight
    
    return nodes[index]; // return the value of the node with the associated index
}


/*///////////////////////
//// Model Parameter ////
////    Structure    ////
///////////////////////*/

struct Parameters
{

    // variables for past model structure
    int past_model = 0; // switch indicator 0: fixed N; 1: poisson N; 2: gamma-poisson N; -- 0
    int N = 0; // fixed past spillovers -- 1
    double lambda = 1.0; // constant rate of spillover -- 2
    double k = 1.0; // gamma shape parameter for spillover rate -- 3
    double theta = 1.0; // gamma scale parameter for spillover rate -- 4
    
    
    // variable for future model structure
    int future_model = 0; // switch indicator 0: fixed M; 1: poisson M; 2: gamma-poisson M; -- 5
    int M = 0; // fixed future spillovers -- 6
    double lambda_f = 1.0; // constant rate of spillover -- 7
    double k_f = 1.0; // gamma shape parameter for spillover rate -- 8
    double theta_f = 1.0; // gamma scale parameter for spillover rate -- 9
    
    
    // model output and other control parameters
    int redraw = 0; // indicator to redraw lambda in the future (0 - no; 1 - yes) -- 10
    int verbose = 0; // verbose statement indicator for main and model function -- 11
    int std_out = 0; // output type {note that writing to file extends model runtime} (0 - file; 1 - console) -- 12
    int output_mode = 0; // should full output be provided or only information necessary for standard plots {effect on runtime has not been fully assessed} (0 - full; 1 - minimal) -- 13
    std::string batchname = "simData_test"; // output file name tag -- 14
    
    // Model prior parameters
    int prior_type = 0; // prior parameter options (0 - beta; 1 - discrete; 2 - B-splines) -- 15
    double a_prior = 1.0; // shape1 beta prior parameter -- 16
    double b_prior = 1.0; // shape2 beta prior parameter -- 17
    int H_crit = 0; // permitted number of past host jumps, ideally 0 -- 18
    
    int tol = 1e9; // simulation tolerance. a poorly informed prior can cause excessive runtime
    
    // special parameters -- read from file created in the R wrapper function
    std::vector<double> nodes; // nodes for discrete prior
    std::vector<double> weights; // corresponding weights for node values. not necessarily normalized -- also used for beta mixture models
    
    std::vector<double> a_vec; // nodes for discrete prior
    std::vector<double> b_vec; // corresponding weights for node values. not necessarily normalized
    
    std::vector<std::string> conf_v;

    //Constructor either reads parameters from standard input (if no argument is given),
    //or from file (argument 1(.
    Parameters(int argc, char* argv[])
    {
        conf_v.reserve(200);
        std::stringstream buffer;
        if(argc == 3) //Config file
        {
            std::ifstream f(argv[2]);
            if(f.is_open())
            {
                buffer << f.rdbuf();
            }
            else
            {
                std::cout << "Failed to read config file \"" << argv[2] << "\"" << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        else
        {
            for(int i=2; i<argc; ++i)
            {
                buffer << argv[i] << std::endl;
            }
        }
        while(!buffer.eof())
        { // until end of the stream
            std::string line = "";
            std::stringstream line_ss;
            // First get a whole line from the config file
            std::getline(buffer, line);
            // Put that in a stringstream (required for getline) and use getline again
            // to only read up until the comment character (delimiter).
            line_ss << line;
            std::getline(line_ss, line, '#');
            // If there is whitespace between the value of interest and the comment character
            // this will be read as part of the value. Therefore strip it of whitespace first.
            line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());
            if(line.size() != 0)
            {
                conf_v.push_back(line);
            }
        }

        size_t conf_length = 19; //number of model input parameters
        if (conf_v.size()!=conf_length)
        {
            std::cout << "Expected configuration file with " << conf_length << " options, loaded file with "
                      << conf_v.size() << " lines." << std::endl;
            exit(EXIT_FAILURE);
        }

        
        // variables for past model structure
        past_model = std::stoi(conf_v[0]); // switch indicator 0: fixed N; 1: poisson N; 2: gamma-poisson N; -- 0
        N = std::stoi(conf_v[1]); // fixed past spillovers -- 1
        lambda = std::stod(conf_v[2]); // constant rate of spillover -- 2
        k = std::stod(conf_v[3]); // gamma shape parameter for spillover rate -- 3
        theta = std::stod(conf_v[4]); // gamma scale parameter for spillover rate -- 4
        
        
        // variable for future model structure
        future_model = std::stoi(conf_v[5]); // switch indicator 0: fixed M; 1: poisson M; 2: gamma-poisson M; -- 5
        M = std::stoi(conf_v[6]); // fixed future spillovers -- 6
        lambda_f = std::stod(conf_v[7]); // constant rate of spillover -- 7
        k_f = std::stod(conf_v[8]); // gamma shape parameter for spillover rate -- 8
        theta_f = std::stod(conf_v[9]); // gamma scale parameter for spillover rate -- 9
        
        
        // model output and other control parameters
        redraw = std::stoi(conf_v[10]); // indicator to redraw lambda in the future (0 - no; 1 - yes) -- 10
        verbose = std::stoi(conf_v[11]); // verbose statement indicator for main and model function -- 11
        std_out = std::stoi(conf_v[12]); // output type {note that writing to file extends model runtime} (0 - file; 1 - console) -- 12
        output_mode = std::stoi(conf_v[13]); // should full output be provided or only information necessary for standard plots {effect on runtime has not been fully assessed} (0 - full; 1 - minimal) -- 13
        batchname = batchname = conf_v[14]; // output file name tag -- 14
        
        // Model prior parameters
        prior_type = std::stoi(conf_v[15]); // prior parameter options (0 - beta; 1 - discrete; 2 - B-splines) -- 15
        a_prior = std::stod(conf_v[16]); // shape1 beta prior parameter -- 16
        b_prior = std::stod(conf_v[17]); // shape2 beta prior parameter -- 17
        H_crit = std::stod(conf_v[18]); // permitted number of past host jumps, ideally 0 -- 18
        
        if(prior_type == 1){ // read input from discrete prior parameter file
            
            // path will need to be corrected for other locations
            std::ifstream parms("./simulation_files/discreteParms.csv");
            // std::ifstream parms("/Users/brandonsimony/Desktop/Repos/HostJump_Model/src/discreteParms.csv");
            if(parms.is_open()){
                
                std::string line = "";
                while(std::getline(parms, line)){
                    
                    double node, weight;
                    std::string tmp = "";
                    std::stringstream iss(line);
                    
                    std::getline(iss, tmp, ',');
                    node = std::stod(tmp); // convert string to double
                    nodes.push_back(node); // add node value to vector
                    tmp = "";
                    
                    std::getline(iss, tmp, ',');
                    weight = std::stod(tmp); // convert string to double
                    weights.push_back(weight); // add weight value to vector
                     
                }

            }else{
                std::cout << "Failed to open file: discreteParms.csv" << std::endl;
                exit(EXIT_FAILURE);
            }
            parms.close();
                
        }else if(prior_type == 2){ // read input parameters from the beta mixture parameter file
            
            // path will need to be corrected for other locations
            std::ifstream parms("./simulation_files/betaMix_pars.csv");
            if(parms.is_open()){
                
                std::string line = "";
                while(std::getline(parms, line)){
                    
                    double a, b, norm_weight, weight;
                    std::string tmp = "";
                    std::stringstream iss(line);
                    
                    std::getline(iss, tmp, ',');
                    a = std::stod(tmp); // convert string to double
                    a_vec.push_back(a); // add node value to vector
                    tmp = "";
                    
                    std::getline(iss, tmp, ',');
                    b = std::stod(tmp); // convert string to double
                    b_vec.push_back(b); // add node value to vector
                    tmp = "";
                    
                    // necessary update to code for poissonExact.cpp script which reads in the same parameter file, but requires normalized weights instead of cumulative weights for computation. Here these normalized weights are unnecessary but the line must be read to ensure proper values are used. 
                    std::getline(iss, tmp, ',');
                    // norm_weight = std::stod(tmp); // convert string to double
                    // nodes.push_back(norm_weight); // add normalized weight value to vector
                    
                    std::getline(iss, tmp, ',');
                    weight = std::stod(tmp); // convert string to double
                    weights.push_back(weight); // add weight value to vector
                     
                }

            }else{
                std::cout << "Failed to read file: betaMix_pars.csv" << std::endl;
                exit(EXIT_FAILURE);
            }
            parms.close();
            
        } // end reading in prior parameter files
        
        
    }
};


// structure to hold "accepted" simulated data sets
struct simData{
    
    std::vector<int> N; // saved number of past spillovers
    std::vector<int> M; // drawn number of future spillovers
    std::vector<int> H_f; // count of future spillovers that resulted in host jumps
    std::vector<double> lambda; // saved drawn or fixed value for lambda
    std::vector<double> phi; // saved drawn value of phi
    std::vector<int> attempt; // indicator for number of trial datasets before a draw with no host jumps in the past
    double total_prob; // fraction of replicates that had a future host jump given no past host jump
    
    simData(int n_reps){
        
        N.assign(n_reps, -1);
        M.assign(n_reps, -1);
        H_f.assign(n_reps, -1);
        lambda.assign(n_reps, -1.0);
        phi.assign(n_reps, -1.0);
        attempt.assign(n_reps, -1);
        total_prob = 0;
        
    }
    
};


/*/////////////////////
// Model Simulations //
/////////////////////*/

simData pastSim( std::ofstream &simResults, std::ofstream &errOut, const int n_reps, const Parameters& p){
    
    simData data(n_reps);
    
    // generate hypothetical past data sets
    int idx = 0;
    while( idx < n_reps ){
        
        int N_draw = p.N; // past spillovers
        int M_draw = p.M; // future spillovers
        double phi_draw = 0; // host jump probability given spillover
        int H_past = 0; // past host jumps
        int H_f = 0; // future host jumps
        double lambda = p.lambda; // spillover rate
        double break_val = 0.0; // alternative loop break condition if simulating exact data is impossible, i.e. N_draw < H_crit always.
        data.attempt[idx] ++;
        if(data.attempt[idx] > p.tol & p.tol > 0){
            std::cout << "No successful datasets simulated in " << p.tol << " attempts. Prior distribution is likely unrealistic for this scenario" << std::endl;
            errOut << "No successful datasets simulated in " << p.tol << " attempts. Prior distribution is likely unrealistic for this scenario" << std::endl;
            errOut.flush();
            exit(EXIT_FAILURE);
        }
        
        // determine appropriate past parameters
        switch(p.past_model){
            
            // simulate past with fixed number of spillover events
            case 0:
                break_val = p.N;
                if(p.N < p.H_crit){break_val = 0.0;} // force alternative break since data simulation will be impossible
                // draw phi from the desired prior distribution
                switch(p.prior_type){
                    // beta prior
                    case 0:
                        phi_draw = draw_beta_gsl(p.a_prior, p.b_prior);
                    break;
                        
                    // discrete prior
                    case 1:
                        phi_draw = draw_discrete_gsl(p.nodes, p.weights);
                    break;
                        
                    // beta mixture prior
                    case 2:
                        phi_draw = draw_beta_mix(p.a_vec, p.b_vec, p.weights);
                    break;
                        
                    default:
                        std::cout << "invalid prior choice: " << p.prior_type << std::endl;
                        exit(EXIT_FAILURE);
                }
                H_past = draw_binom_gsl(N_draw, phi_draw);
                
            break;
                
            // simulate past with poisson number of spillovers -- using default lambda
            case 1:
                
                break_val = p.lambda; // force alternative break if lambda == 0
                // draw phi from the desired prior distribution
                switch(p.prior_type){
                    // beta prior
                    case 0:
                        phi_draw = draw_beta_gsl(p.a_prior, p.b_prior);
                    break;
                        
                    // discrete prior
                    case 1:
                        phi_draw = draw_discrete_gsl(p.nodes, p.weights);
                    break;
                        
                    // beta mixture prior
                    case 2:
                        phi_draw = draw_beta_mix(p.a_vec, p.b_vec, p.weights);
                    break;
                        
                    default:
                        std::cout << "invalid prior choice: " << p.prior_type << std::endl;
                        exit(EXIT_FAILURE);
                }
                N_draw = draw_poisson_gsl(lambda);
                H_past = draw_binom_gsl(N_draw, phi_draw);
                
            break;
                
            // simulate past with gamma-poisson number of spillovers -- resets lambda default
            case 2:
                
                break_val = std::min(p.k, p.theta); // force alternative break if either parameter is 0
                // draw phi from the desired prior distribution
                switch(p.prior_type){
                    // beta prior
                    case 0:
                        phi_draw = draw_beta_gsl(p.a_prior, p.b_prior);
                    break;
                        
                    // discrete prior
                    case 1:
                        phi_draw = draw_discrete_gsl(p.nodes, p.weights);
                    break;
                        
                    // beta mixture prior
                    case 2:
                        phi_draw = draw_beta_mix(p.a_vec, p.b_vec, p.weights);
                    break;
                        
                    default:
                        std::cout << "invalid prior choice: " << p.prior_type << std::endl;
                        exit(EXIT_FAILURE);
                }
                lambda = draw_gamma_gsl(p.k, p.theta);
                N_draw = draw_poisson_gsl(lambda);
                H_past = draw_binom_gsl(N_draw, phi_draw);
                
            break;
                
            default:
                std::cout << "invalid case choice: " << p.past_model << std::endl;
                exit(EXIT_FAILURE);
                
        }
        
        
        // simulate future conditional on some number of past host jumps
        if(H_past == p.H_crit | break_val <= 0){
            
            // save values that satisfy no past host jumps
            data.N[idx] = N_draw;
            data.lambda[idx] = lambda;
            data.phi[idx] = phi_draw;
            
            // determine appropriate past parameters
            switch(p.future_model){
                
                // simulate past with fixed number of spillover events
                case 0:
                    H_f = draw_binom_gsl(M_draw, phi_draw);
                break;
                    
                // simulate past with poisson number of spillovers --
                case 1:
                    if(p.redraw == 1 || lambda == 0.0){lambda = p.lambda_f;} // uses future value of lambda if past lambda was undefined or 0 or if redraws are requested
                    M_draw = draw_poisson_gsl(lambda);
                    H_f = draw_binom_gsl(M_draw, phi_draw);
                break;
                    
                // simulate past with gamma-poisson number of spillovers
                case 2:
                    if(p.redraw == 1 || lambda == 0.0){lambda = draw_gamma_gsl(p.k_f, p.theta_f);} // updates future value of lambda if past lambda was undefined or 0 or if redraws are requested
                    M_draw = draw_poisson_gsl(lambda);
                    H_f = draw_binom_gsl(M_draw, phi_draw);
                break;
                    
                default:
                    std::cout << "invalid case choice: " << p.future_model << std::endl;
                    exit(EXIT_FAILURE);
                    
            }
            
            if(H_f > 0){ data.total_prob = data.total_prob + 1; }
            data.M[idx] = M_draw;
            data.H_f[idx] = H_f;
            idx ++;
        }
        
    } // end while loop
    
    data.total_prob = double( double(data.total_prob) / double(n_reps) );
    return(data);
    
}

/*////////////
//// Main ////
////////////*/

int main(int argc, char* argv[]){
    
    /*
     If you run this from the command line, you should provide two arguments - the
     number of replicates and a config file (argc==3). In that case the results will
     be output to a file. The alternative is to provide all the parameters in a single
     long string instead of the config file. This is done by the R-wrapper. In that
     case the results will be written to a .txt file and saved to a table in R.
     */
    
   
    if(argc < 3)
    {
        std::cout << "Usage: " << argv[0] << " <n replicates> <config file / string of parameter values>" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    // initialize parameters read from command line input
    int n_reps = std::stoi(argv[1]);
    Parameters p(argc, argv);
    
    
    // open full output file
    std::string fname = "./src/simulation_files/" + p.batchname + "_fullSim.txt";
    std::ofstream simResults(fname);
    //assert(simResults.is_open());
    
    // open error file
    std::string err_fname = "./src/simulation_files/" + p.batchname + "_err.txt";
    std::ofstream errOut(err_fname);
    // assert(errOut.is_open());
    
    // end open output files
    
    
    simData data = pastSim( simResults, errOut, n_reps, p);
    
    // model output blocks -- either write to file or console output
    if(p.std_out == 0){
        
        switch(p.output_mode){
                
            // full output
            case 0:
                
                // write column labels
                simResults << "total_prob;" << "lambda;" << "phi;" << "N;" << "M;" << "H_future" << std::endl;
                simResults.flush();
                
                for(int i = 0; i < n_reps; i++){
                    // write results from each individual simulation
                    simResults << data.total_prob << ";" << data.lambda[i] << ";" << data.phi[i] << ";" << data.N[i] << ";" << data.M[i] << ";" << data.H_f[i] << std::endl;
                    simResults.flush();
                }
                
            break;
                
            // minimal output
            case 1:
                simResults << "total_prob" << std::endl;
                simResults << data.total_prob << std::endl;
                simResults.flush();
            break;
                
            // minimal output also as default
            default:
                simResults << "total_prob" << std::endl;
                simResults << data.total_prob << std::endl;
                simResults.flush();
        }
        
    }else{
        
        switch(p.output_mode){
                
            // full output
            case 0:
                
                // column labels
                std::cout << "total_prob;" << "lambda;" << "phi;" << "N;" << "M;" << "H_future" << std::endl;
                for(int i = 0; i < n_reps; i++){
                    // print out results from each individual simulation
                    std::cout << data.total_prob << ";" << data.lambda[i] << ";" << data.phi[i] << ";" << data.N[i] << ";" << data.M[i] << ";" << data.H_f[i] << std::endl;
                }
                
            break;
                
            // minimal output
            case 1:
                std::cout << "total_prob" << std::endl;
                std::cout << data.total_prob << std::endl;
            break;
                
            // minimal output also as default
            default:
                std::cout << "total_prob" << std::endl;
                std::cout << data.total_prob << std::endl;
        }
        
    }
          
    simResults.close(); // result file
    errOut.close(); // error file
    
    return 0;
    
}


