#include "swap.h"

double dXCM, dYCM;
int dataCounter=0;
int cycle=0;
// double swapCount[N] = {0};

// Monte Carlo Simulation
void MC(std::string out, int n_log, int n_lin){
    int cycleCounter = 0;
    double deltaX[N], deltaY[N], deltaR2[N], R2Max = 0;
    // Building snapshots list (log-spaced)
    std::vector < std::pair <double, double>> pairs;
    std::vector <double> samplePoints, twPoints;
    double endingPoints[cycles], linPoints[n_lin];
    double exponents = log10(tau)/(n_log-1);
    std::string in = out + "configs/";
    for(int c=0; c<cycles; c++){
        for (int x = 0; x < n_log; x++){
            double value = tw*c + floor(pow(10,exponents*(x)));
            std::pair <double,double> p = {value, c};
            int f = std::count(pairs.begin(), pairs.end(), p);
            if(f==0){
                pairs.emplace_back(value, c);
            // this if condition is actually relevent because of the floor function
            }
        }
    }

    // Sorting
    std::sort(pairs.begin(), pairs.end());
    for (auto p: pairs){
        samplePoints.push_back(p.first); twPoints.push_back(p.second);
    }

    // Ending points
    for(int c=0;c<cycles;c++){
        endingPoints[c] = c*tw + tau;
    }
    // Linspaced points
    for (int k=1;k<=n_lin;k++){
        linPoints[k] = (tau/(n_lin))*k;
    }
    // File writing
    std::ofstream log_obs, log_cfg; 
    // log_ploc, log_p;
    // log_sigma;
    std::string out_cfg = out + "configs/";
    // std::string out_ploc = out + "micro_corr/";
    // std::string out_sigma = out + "sigma_scan/";
    log_obs.open(out + "obs.txt");
    log_obs << "t" << " " << "cycle";
    for (std::string obs: allObs){
        log_obs << " " << obs;
    } log_obs << std::endl;
    // log_p.open(out + "space_corr.txt");
    log_obs << std::scientific << std::setprecision(8);
    // log_p << std::scientific << std::setprecision(8);
    // creating outdir if not existing
    // fs::create_directory(out_ploc);

    for(double t_: samplePoints){
        // auto start = std::chrono::high_resolution_clock::now();
        int t = int(t_);
        std::string cfg = in + "cfg_" + std::to_string(t) + ".xy";
        ReadCFG(cfg);
        
        if(t==1){
            for(int i=0;i<N;i++){
                Xref[i] = Xfull[i]; Yref[i] = Yfull[i];
            }
        }

        // Updating reference observables
        if((t-1)%tw == 0 && cycleCounter < cycles){
            UpdateAge(cycleCounter); cycleCounter++;
        } 
        
        UpdateNL(); UpdateNN(); // UpdateRL(); // updating nearest neighbours
        
        dXCM = 0; dYCM = 0;
        if(t!=1){
            for (int i=0;i<N;i++){
                double dX = Xfull[i]-Xref[i], dY = Yfull[i]-Yref[i];
                dXCM += dX; dYCM += dY;
            } dXCM /= N; dYCM /= N;
        }
        // int cycle = twPoints[dataCounter];
        // Configs
        // log_ploc.open(out_ploc + "products_loc_" + std::to_string(t) + ".txt");
        // log_ploc << std::scientific << std::setprecision(8);
        // for (int i = 0; i<N; i++){
        //     std::vector <double> disp_loc = MicroDispCorrLoc(i);
        //     // std::vector <double> u_sigma = SigmaScan(i);
        //     for (int k=0;k<nr;k++){
        //         log_ploc << disp_loc[k] << " ";
        //     } log_ploc << std::endl;
        // };
        // log_ploc.close();

        log_obs << t << " " << cycle;
                for (std::string obs: allObs){
                    log_obs << " " << whichObs(obs, cycle);
                } log_obs << std::endl;

        // log_p << t << " ";
        // std::vector <double> disp = MicroDispCorr();
        // for (int k=0;k<nr;k++){
        //         log_p << disp[k] << " ";
        //     } log_p << std::endl;
        // saving format: timestep Vtot MSD Fs CB 

        // dataCounter++;

        // auto end = std::chrono::high_resolution_clock::now();
        // std::chrono::duration<double> duration = end - start;
        // std::cout << duration.count() << std::endl;
        std::cout << t << std::endl;
    }
    log_obs.close();
    // log_p.close();
}

//  Tries displacing one particle j by vector dr = (dx, dy)
void TryDisp(int j){
    double dx = (ranf()-0.5)*deltaMax;
    double dy = (ranf()-0.5)*deltaMax;
    double Xnew = Pshift(X[j]+dx);
    double Ynew = Pshift(Y[j]+dy);
    double deltaE = V(Xnew, Ynew, S[j], j) - V(X[j], Y[j], S[j], j);
    // why is the modulus function not in deltaE ?
    if (deltaE < 0){
        // Xnew = fmod(X[j],Size);
        X[j] = Xnew; //Check modulus function
        Y[j] = Ynew;
        Xfull[j] = Xfull[j]+dx;
        Yfull[j] = Yfull[j]+dy;
    }
    else if (exp(-deltaE/T) > ranf()){
        X[j] = Xnew;
        Y[j] = Ynew;
        Xfull[j] = Xfull[j]+dx;
        Yfull[j] = Yfull[j]+dy;
    }
}

//  Tries swapping the pair of particles j, k
void TrySwap(int j, int k){
    double deltaS = std::abs (S[j]-S[k]);
    if(deltaS<=deltaSMax){
        double deltaE = V(X[j],Y[j],S[k],j)+V(X[k],Y[k],S[j],k)-V(X[j],Y[j],S[j],j)-V(X[k],Y[k],S[k],k);
        if (deltaE < 0){
            double Rnew = S[k];
            S[k] = S[j];
            S[j] = Rnew;
            // swapCount[j] += 1; swapCount[k] += 1;
        }
        else if (exp(-deltaE/T) > ranf()){
            double Rnew = S[k];
            S[k] = S[j];
            S[j] = Rnew;
            // swapCount[j] += 1; swapCount[k] += 1;
        }
    } else{
        // pass
    }
}

double whichObs(std::string obs, int cycl){
    if (obs=="MSD") return MSD();
    else if (obs=="U") return VTotal()/(2*N);
    else if (obs=="Cb") return CB(cycl);
    else if (obs=="Fs") return FS(cycl);
}