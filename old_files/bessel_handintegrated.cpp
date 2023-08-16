#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include "eigen/Eigen/Dense"
#include "variables.h"
#include <chrono>
#include <fstream>
#include <string>
#include <sstream>


const double hbar = 1; /*1.0545718e-34;*/
const double pi = 3.14159265358979323846;
const double sqrtPi = 1.7724538509055159;

const double length = 1; /*długość małego walca == L*/
const double radius = 1; /*promień małego walca == sigma*/
const double mass = 1;

/*przy 2e3 liczenie trwa godzinę dla 515 elementów, 7 sekund na każdy*/
const int NumberOfSteps = 2e3;
const double dstep=1./double(NumberOfSteps);


double bessel_polynomial(double rho, int l, int m){
    return std::cyl_bessel_j(double(l),double(bessel_zeros.at(l-1).at(m)*rho/radius));
}

double delta_integral(int l1, int l2, int m){
    double integrated_value = 0;
    for (int i0=0 ; i0<NumberOfSteps ; i0++){
        double rho=i0*dstep*radius;
        integrated_value +=
            bessel_polynomial(rho, l1, m)*bessel_polynomial(rho, l2, m)*dstep*radius;
    }
    return integrated_value;
}

double laplace_integral(int l1, int l2, int m){
    double integrated_value = 0;
    for (int i0=0 ; i0<NumberOfSteps ; i0++){
        double rho1=i0*dstep*radius+1e-8;
        for (int i1=0 ; i1<NumberOfSteps ; i1++){
            double rho2=i1*dstep*radius+1e-8;
            integrated_value +=
            bessel_polynomial(rho1, l1, m)*bessel_polynomial(rho2, l2, m)*(rho1*rho1+rho2*rho2)/(rho1*rho1*rho2*rho2)*dstep*radius*dstep*radius;
        }
    }
    return integrated_value;
}

double combined_integral(double g, int n1, int l1, int n2, int l2, int m){
    double dlta_integral_value = 0;
    double lplc_integral_value = 0;
    for (int i0=0 ; i0<NumberOfSteps ; i0++){
        double rho1=i0*radius+1e-8;
        dlta_integral_value +=
            bessel_polynomial(rho1, l1, m)*bessel_polynomial(rho1, l2, m);
        for (int i1=0 ; i1<NumberOfSteps ; i1++){
            double rho2=i1*radius+1e-8;
            lplc_integral_value +=
                bessel_polynomial(rho1, l1, m)*bessel_polynomial(rho2, l2, m)*(rho1*rho1+rho2*rho2)/(rho1*rho1*rho2*rho2);
        }
    }
    lplc_integral_value*=dstep*radius*dstep*radius;
    dlta_integral_value*=dstep*radius;

    return -8*pi*pi*(n1*n1+n2*n2)*lplc_integral_value+2*pi*g*length*dlta_integral_value;
}


int main(int argc, char** argv){
auto start = std::chrono::high_resolution_clock::now();

    double g=1.;

    std::cout << "number of steps: " << NumberOfSteps << std::endl;

/*
    std::cout <<
    8*pi*pi*(n1*n1+n2*n2)*laplace_integral(l1, l2, m)+2*pi*g*length*delta_integral(l1, l2, m)
    <<std::endl;
*/


    std::vector<std::string> parameters = {"n1", "l1", "n2", "l2", "m", "integral_value"};
    std::vector<std::vector<double>> values;

    for (int n1 = 0; n1 <= 1; ++n1) {
        for (int n2 = 0; n2 <= 1; ++n2) {
            for (int l1 = 1; l1 <= 1; ++l1) {
                for (int l2 = 1; l2 <= 1; ++l2) {
                    for (int m = 0; m <= 1; ++m) {
                        std::cout << "n1: " << n1 << ", n2: " << n2 << ", l1: " << l1 
                                  << ", l2: " << l2 << ", m: " << m << ", element tensora: ";

                        
                        double tensor_element = combined_integral(g, n1, l1, n2, l2, m);
                        /*
                        double tensor_element = 8*pi*pi*(n1*n1+n2*n2)*laplace_integral(l1, l2, m)+2*pi*g*length*delta_integral(l1, l2, m);
                        */
                        std::cout << tensor_element << std::endl;

                        values.push_back(std::vector<double>({double(n1),double(n2),double(l1),double(l2),double(m), tensor_element}));

                    }
                }
            }
        }
    }

auto end = std::chrono::high_resolution_clock::now();
double elapsed_seconds = std::chrono::duration<double>(end - start).count();
std::cout << "Elapsed time: " << elapsed_seconds << " seconds." << std::endl;


std::ostringstream oss;
oss << "Elapsed time: " << elapsed_seconds << "\tg = " << g;
std::string extrainfo = oss.str();
    functions Function;
    Function.write_to_file("output.txt", parameters, values, extrainfo);

    return 0;
}

/*
działa to strasznie dziwnie, nawet nie pozwala na symetryzację. Może popełniłem błąd w funkcji A?
*/