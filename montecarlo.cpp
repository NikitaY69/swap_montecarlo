#include "swap.h"

// Monte Carlo Simulation
void MC(std::string out, int ss){
    int dataCounter = 0, cycleCounter = 0;
    double deltaX[N], deltaY[N], deltaR2[N], R2Max = 0;

    // Building snapshots list (log-spaced)
    std::vector < std::pair <double, double>> pairs;
    std::vector <double> samplePoints, twPoints;
    double exponents = log10(tau)/(ss-1);

    for(int c=0; c<cycles; c++){
        for (int x = 0; x < ss; x++){
            double value = tw*c + floor(pow(10,exponents*(x)));
            std::pair <double,double> p = {value, c};
            int f = std::count(pairs.begin(), pairs.end(), p);
            // int f = Find(samplePoints, value);
            // std::count(samplePoints.begin(), samplePoints.end(), value);
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

    // File writing
    std::ofstream log_obs, log_cfg;
    log_obs.open(out + "obs.txt");
    log_cfg.open(out + "cfg.xy");
    log_obs << std::scientific << std::setprecision(8);
    log_cfg << std::scientific << std::setprecision(8);

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
                NL.clear();
                UpdateList();
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
        log_cfg << S[451] << " " << X[233] << " " << Y[1432] << std::endl;      
        log_obs << VTotal()/(2*N) << std::endl;          
        
        // Doing the MC
        for (int i = 0; i < N; i++){
            if (ranf() > 0.2) TryDisp(i); //Displacement probability 0.8
            else TrySwap(i,floor(ranf()*N)); //Swap probability 0.2
        }

        if((t-1)%100==0) std::cout << (t-1) << std::endl;; // Counting steps
    };
    log_obs.close(), log_cfg.close();
}

//  Tries displacing one particle j by vector dr = (dx, dy)
void TryDisp(int j){
    double dx = (ranf()-0.5)*deltaMax;
    double dy = (ranf()-0.5)*deltaMax;
    double deltaE = V(X[j] + dx, Y[j] + dy, S[j], j) - V(X[j], Y[j], S[j], j);
    // why is the modulus function not in deltaE ?
    if (deltaE < 0){
        X[j] = fmod((X[j]+dx),Size); //Check modulus function
        Y[j] = fmod((Y[j]+dy),Size);
        Xfull[j] = Xfull[j]+dx;
        Yfull[j] = Yfull[j]+dy;
    }
    else if (exp(-deltaE/T) > ranf()){
        X[j] = fmod((X[j]+dx),Size);
        Y[j] = fmod((Y[j]+dy),Size);
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