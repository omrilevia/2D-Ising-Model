#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <random>
#include <chrono>
#include <unistd.h>
#include <time.h>
using namespace std;

ofstream Ising_data("Ising_data.txt",ios::out);

const int L = 32;
double J = 1.0;
int lattice[L][L];
int delta_m = 0;
double T = 5.0;
double dT = 0.1;
double T_min = 0.5;

int int_samples = 10000;

double N_and_site = 1.0/(double(int_samples*L*L));

/*  *  *  *  *  *  *  *  *  *  *  *
 
 * Different type of random number generators to use
 
 */

unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator (seed);
std::uniform_real_distribution<double> dist_1(0.0,1.0);
std::uniform_int_distribution<int> L_dist(0,L-1);


//initialize lattice

void init_lattice(int lattice[L][L]){
    for (int i = 0; i < (L); i++){
        for (int j = 0; j < (L); j++){
            
            double r = dist_1(generator);
            if (r<=0.5){
                lattice[i][j] = 1;
            }else{
                lattice[i][j] = -1;
            }
         //cout<<lattice[i][j]<<endl;
        }
    }

}

//main algorithm
int ising_MC(int lattice[L][L],double T,int& delta_m){
    //find a random site
    int i = L_dist(generator);
    usleep(15);
    int j = L_dist(generator);
    //cout<<i<<"  "<<j<<endl;
    //implement periodic boundary conditions, compute energy on that site
    int up,down,left,right,energy;
    up = (j == (L-1)) ? 0 : j+1;
    down = (j==0) ? L-1 : j-1;
    left = (i==0) ? L-1 : i-1;
    right = (i==(L-1)) ? 0 : i+1;
    
    energy=-1*(lattice[i][j])*(lattice[left][j]+lattice[right][j]+lattice[i][up]+lattice[i][down]);
    //cout<<energy<<endl;
    //delta_E from flipping the i,j spin is:
    int delta_e = -2*energy;
    double r = dist_1(generator);
    
    if(delta_e<0){
        lattice[i][j] = -1*lattice[i][j];
        delta_m = 2*lattice[i][j];
    }else if(r<=exp(-1.0*double(delta_e/T))){
        lattice[i][j] = -1*lattice[i][j];
        delta_m = 2*lattice[i][j];
    }else{
        delta_e = 0;
        delta_m = 0;
    
    }return delta_e;

}

//function for alleviation of transient effects
void transient_effects(int& delta_m){
    for (int i=0;i<1000;i++){
        for(int j=0;j<100;j++){
            ising_MC(lattice,T,delta_m);
        }
    }
}

//function for calculating energy at a single site
int get_energy(int lattice[L][L],int i,int j){
    
    //implement periodic boundary conditions
    int up,down,left,right,energy;
    up = (j == (L-1)) ? 0 : j+1;
    down = (j==0) ? L-1 : j-1;
    left = (i==0) ? L-1 : i-1;
    right = (i==(L-1)) ? 0 : i+1;
    
    energy=-1*(lattice[i][j])*(lattice[left][j]+lattice[right][j]+lattice[i][up]+lattice[i][down]);
    
    return energy;
    
}

//function for computing total energy of lattice
int energy_total(){
    int energy = 0;
    for(int i =0;i<L;i++){
        for(int j=0;j<L;j++){
            energy += get_energy(lattice,i,j);
        }
    }
    return energy;
}

//function for computing total magnetization of lattice
int magnet_total(){
    int m=0;
    for (int i=0;i<(L);i++){
        for(int j=0;j<(L);j++){
            m+=lattice[i][j];
        }
    }
    return m;
}

int main(){
    //declaring variables to be used in calculating the observables
    double E=0.0,E_sqr=0.0,E_sqr_avg=0.0,E_avg=0.0,e_tot=0.0,e_tot_sqr=0.0;
    double M=0.0,M_sqr=0.0,M_sqr_avg=0.0,M_avg=0.0,m_tot=0.0,m_tot_sqr=0.0;
    double Mabs=0.0,Mabs_avg=0.0,Mq_avg=0.0,m_abs_tot=0.0;
    
    
    init_lattice(lattice);
    
    //Temperature loop
    do{
        transient_effects(delta_m);
        
        //assign observables
        M = magnet_total();
        Mabs = abs(magnet_total());
        E = energy_total();
//
        //variables to be summed are initialized at beginning of temperature iteration
        e_tot = 0.0;
        e_tot_sqr = 0.0;
        m_tot = 0.0;
        m_tot_sqr = 0.0;
        m_abs_tot = 0.0;
        
        //Monte Carlo Loop
        for (int i=0;i<int_samples;i++){
            
            //Metropolis loop
            for (int j=0;j<100;j++){
                
                E += 2.0*ising_MC(lattice,T,delta_m);
                M += delta_m;
                Mabs += abs(delta_m);
            }
            //cout<<E<<endl;
            //sum observables
            e_tot += E/2.0;
            e_tot_sqr += (E/2.0)*(E/2.0);
            m_tot += M;
            m_tot_sqr += M*M;
            
            m_abs_tot += pow(M*M,0.5);
        }

        //average observables
        E_avg = e_tot*N_and_site;
        E_sqr_avg = e_tot_sqr*N_and_site;
        M_avg = m_tot*N_and_site;
        M_sqr_avg = m_tot_sqr*N_and_site;
        Mabs_avg = m_abs_tot*N_and_site;
        
        
        //output to file
        Ising_data<<T<<                                                //temperature T
        "     "<<M_avg<<"     "<<Mabs_avg<<"     "<<M_sqr_avg<<        //<M>;<|M|>;<M^2> per spin
        "     "<<(M_sqr_avg-(M_avg*M_avg*L*L))/T<<                     //susceptibility per spin X
        "     "<<(M_sqr_avg-(Mabs_avg*Mabs_avg*L*L))/T<<               //susceptibility per spin X
        "     "<<E_avg<<"     "<<E_sqr_avg<<                           //<E>;<E^2> per spin
        "     "<<(E_sqr_avg-(E_avg*E_avg*L*L))/(T*T)<<endl;            //heat capacity C per spin
        

        //decrease T
        T = T-dT;

    

        

    }while(T>=T_min);


    
}
