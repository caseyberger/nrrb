//
// Created by Casey Berger on 2/5/18.
//

#include "test.h"
#include <iostream>
#include <vector>

void print_mu_vals(std::vector<double> vec){
    std::vector<double> mu_vals = vec;
    int n_mu;
    n_mu = (int) mu_vals.size();
    std::cout << "mu values to compute are: "<< std::endl;
    for(int i =0; i<n_mu;i++){
        std::cout << mu_vals[i];
        if(i < n_mu-1){std::cout << ",";}
        else{std::cout << std::endl;}
    }
}