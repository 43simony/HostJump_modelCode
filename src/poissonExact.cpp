//
//  poissonExact.cpp
//  
//
//  Created by Brandon Simony on 12/9/24.
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
#include <gsl/gsl_sf_hyperg.h>


struct Parameters
{

    int K = 1;
    double lambda_p = 1.0;
    double lambda_f = 1.0;

    std::vector<double> a_vec; // vector of a_j prior parameters
    std::vector<double> b_vec; // vector of b_j prior parameters
    std::vector<double> norm_weights; // normalized weight for mixture component j
    
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

        K = std::stoi(conf_v[0]); // read in number of mixture components to define expected parameter vector size
        
        size_t conf_length = 3+3*K; //number of model input parameters
        if (conf_v.size()!=conf_length)
        {
            std::cout << "Expected configuration file with " << conf_length << " options, loaded file with "
                      << conf_v.size() << " lines." << std::endl;
            exit(EXIT_FAILURE);
        }

        // variables for past model structure
        lambda_p = std::stod(conf_v[1]);
        lambda_f = std::stod(conf_v[2]);
        
        // assign beta mixture parameter vectors
        for(int i=0; i<K; i++){
            a_vec.push_back(std::stod(conf_v[i+3]));
            b_vec.push_back(std::stod(conf_v[i+K+3]));
            norm_weights.push_back(std::stod(conf_v[i+2*K+3]));
        }
        
        
        /*
        // path will need to be corrected for other locations
        std::ifstream parms("./simulation_files/betaMix_pars.csv");
        if(parms.is_open()){
            
            std::string line = "";
            while(std::getline(parms, line)){
                
                double a, b, cum_weight, norm_weight;
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
                
                std::getline(iss, tmp, ',');
                cum_weight = std::stod(tmp); // convert string to double
                cum_weights.push_back(cum_weight); // add weight value to vector
                
                std::getline(iss, tmp, ',');
                norm_weight = std::stod(tmp); // convert string to double
                norm_weights.push_back(norm_weight); // add weight value to vector
                 
            }

        }else{
            std::cout << "Failed to read file: betaMix_pars.csv" << std::endl;
            exit(EXIT_FAILURE);
        }
        parms.close();
        */
        
    }
};



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
    double t = std::stod(argv[1]);
    Parameters p(argc, argv);
    
    double Denom = 0.0;
    for(int i=0; i<p.a_vec.size(); i++){
        Denom += p.norm_weights[i]*( gsl_sf_hyperg_1F1(p.a_vec[i], p.a_vec[i] + p.b_vec[i], -1.0*p.lambda_p));
    }
    
    double val = 0.0;
    for(int i=0; i<p.a_vec.size(); i++){
        val += p.norm_weights[i]*( gsl_sf_hyperg_1F1( p.a_vec[i], p.a_vec[i] + p.b_vec[i], -1.0*(p.lambda_f+p.lambda_p)) ) / Denom;
    }
    double prob = 1.0 - val;
    
    std::cout << "prob" << std::endl;
    std::cout << prob << std::endl;
    
    
    
    return 0;
    
}


