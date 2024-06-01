#include "swap.h"
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;

// Run parameters
const int tau = 50000;
const int tw = 50000;
const int cycles = 1;
const int steps = tw*(cycles-1)+tau;
const double T = 1; 
std::string motherdir = "/home/allaglo/benchmarks/";
double rSkin, rNL, RUpdate;

// Snapshots
const int dataPoints = 10;

// Initialization of external variables
double X[N], Y[N], S[N], X0[N], Y0[N];
double Xfull[N], Yfull[N], Xref[N], Yref[N];
std::vector < std::array <double, N>> Xtw, Ytw;
std::vector < std::vector<int> > NL, nn_0;
std::vector < std::vector < std::vector <int>>> nn_tw;
//-----------------------------------------------------------------------------
//  main.cpp
int main(int argc, const char * argv[]) {
    
    // User-defined variables
    srand(45128); //Random number generator
    std::string input = motherdir + argv[1];
    std::string outdir = motherdir + argv[2] + "results/";
    rSkin = atof(argv[3]);
    rNL = pow(rC+rSkin,2);
    RUpdate = pow(rSkin,2)/4;

    std::cout << RUpdate << std::endl;
    fs::path out_path = outdir;
    if(!fs::is_directory(out_path)){
        // creating outdir if not existing
        fs::create_directory(outdir);
    }
    
    // Read init config
    std::string line;
    std::ifstream input_file(input);
    if (input_file.is_open()){
        int i = 0; // particle index
        std::vector<std::vector<double>> cfg; // array of configurations
        while (std::getline(input_file, line)){
            double value;
            std::stringstream ss(line);

            cfg.push_back(std::vector<double>());
            while (ss >> value){
                cfg[i].push_back(value);
            }
            S[i] = cfg[i][0]; X[i] = cfg[i][1]; Y[i] = cfg[i][2];
            X0[i] = X[i]; Xfull[i] = X[i]; Xref[i] = X[i];
            Y0[i] = Y[i]; Yfull[i] = Y[i]; Yref[i] = Y[i];
            i++;}
        input_file.close();

    } else {
        std::cout << input << std::endl;
        return 0;
    }

    // // Building list of first neighbours
    for (int i=0; i<N; i++){
        nn_0.push_back(nearest_neighbours(i, x_max));
    }
    UpdateList();

    // // Do simulation with timer
    // double t0 = time(NULL); // Timer
    MC(outdir, dataPoints); 
    // std::cout << rSkin << " " << (time(NULL) - t0) << std::endl; 
    return 0;
}