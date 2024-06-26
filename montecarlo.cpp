#include "swap.h"

double dXCM, dYCM;

// Monte Carlo Simulation
void MC(std::string out, int ss){
    int dataCounter = 0, cycleCounter = 0;
    double deltaX[N], deltaY[N], deltaR2[N], R2Max = 0;
    // Building snapshots list (log-spaced)
    std::vector < std::pair <double, double>> pairs;
    std::vector <double> samplePoints, twPoints;
    double endingPoints[cycles], linPoints[ss];
    double exponents = log10(tau)/(ss-1);

    for(int c=0; c<cycles; c++){
        for (int x = 0; x < ss; x++){
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
    for (int k=1;k<=ss;k++){
        linPoints[k] = (tau/(ss))*k;
    }
    // File writing
    std::ofstream log_obs, log_cfg, log_ploc, log_p;
    std::string out_cfg = out + "configs/";
    std::string out_ploc = out + "micro_p/";
    log_obs.open(out + "obs.txt"), log_p.open(out + "products.txt");
    log_obs << std::scientific << std::setprecision(8);
    log_p << std::scientific << std::setprecision(8);
    // creating outdir if not existing
    fs::create_directory(out_cfg); fs::create_directory(out_ploc);

    for(int t = 1; t <= steps; t++){

        // Updating NL
        if((t-1) % 150 == 0) {//Change number?
            // every 150 steps we check if we need to update the NL
            for (int i = 0; i < N; i++){
                deltaX[i] = bcs(X[i],X0[i]);
                deltaY[i] = bcs(Y[i],Y0[i]);
                deltaR2[i] = deltaX[i]*deltaX[i] + deltaY[i]*deltaY[i];
            R2Max = std::max_element(deltaR2,deltaR2+N)[0];
            }
            if(R2Max > RUpdate){
                UpdateNL();
                R2Max = 0;
                for(int j = 0; j < N; j++){
                    X0[j] = X[j];
                    Y0[j] = Y[j];
                }
            }
        }
    
        // Updating reference observables
        if((t-1)%tw == 0 && cycleCounter < cycles){
            UpdateAge(cycleCounter); cycleCounter++;
        } 

        // Writing observables to text file
        int f = std::count(linPoints, linPoints+ss, 1.0*t);
        // int f = std::count(samplePoints.begin(), samplePoints.end(), t*1.0);
        if(f>0){ // checking if saving time
            UpdateNN(); UpdateRL(); // updating nearest neighbours
            dXCM = 0; dYCM = 0;
            for (int i=0;i<N;i++){
                double deltaX = Xfull[i]-Xref[i], deltaY = Yfull[i]-Yref[i];
                dXCM += deltaX; dYCM += deltaY;
            } dXCM /= N; dYCM /= N;
            for(int s=0; s<f; s++){
                // looping different eventual tws
                int cycle = twPoints[dataCounter];
                if(cycles == 1){
                    // Configs
                    log_cfg.open(out_cfg + "cfg_" + std::to_string(t) + ".xy");
                    log_ploc.open(out_ploc + "products_loc_" + std::to_string(t) + ".txt");
                    log_cfg << std::scientific << std::setprecision(8);
                    log_ploc << std::scientific << std::setprecision(8);
                    for (int i = 0; i<N; i++){
                        std::vector <double> disp_loc = MicroDispCorrLoc(i);
                        log_cfg << S[i] << " " << Xfull[i] << " " << Yfull[i] << std::endl;
                        for (int k=0;k<nr;k++){
                            log_ploc << disp_loc[k] << " ";
                        } log_ploc << DispCorrLoc(i) << std::endl;
                    }
                    log_cfg.close(), log_ploc.close();
                    log_obs << t << " " << MSD() << " " << DispCorr() << " "<< std::endl;
                    log_p << t << " ";
                    std::vector <double> disp = MicroDispCorr();
                    for (int k=0;k<nr;k++){
                            log_p << disp[k] << " ";
                        } log_p << std::endl;
                    // saving format: timestep Vtot MSD Fs CB 
                    
                } else{
                    log_obs << t << " " << cycle << " " << VTotal()/(2*N) << " " <<
                    FS(cycle) << " " << CB(cycle) << std::endl;
                    // saving format: timestep Fs CB 
                }
                dataCounter++;
            }  
        }
        if (cycles > 1 && std::count(endingPoints, endingPoints+cycles, 1.0*t) > 0){
            log_cfg.open(out + "cfg_" + std::to_string(t) + ".xy");
            log_cfg << std::scientific << std::setprecision(8);
            for (int i = 0; i<N; i++){
                log_cfg << S[i] << " " << Xfull[i] << " " << Yfull[i] << std::endl;
            }
            log_cfg.close();
        }
        // Doing the MC
        for (int i = 0; i < N; i++){
            if (ranf() > 0.2) TryDisp(i); //Displacement probability 0.8
            else TrySwap(i,floor(ranf()*N)); //Swap probability 0.2
        }

        if((t-1)%100==0) std::cout << (t-1) << std::endl;; // Counting steps
    };
    log_obs.close(), log_p.close();
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
        }
        else if (exp(-deltaE/T) > ranf()){
            double Rnew = S[k];
            S[k] = S[j];
            S[j] = Rnew;
        }
    } else{
        // pass
    }
}