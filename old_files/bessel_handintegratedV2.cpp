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

typedef long double Double;

const long double hbar = 1e-10; /*1.0545718e-34;*/
const double pi = 3.14159265358979323846;
const double sqrtPi = 1.7724538509055159;

const double length = 1; /*długość małego walca == L*/
const double radius = 1; /*promień małego walca == sigma*/
const double mass = 10;

/*przy 2e3 liczenie trwa godzinę dla 515 elementów, 7 sekund na każdy*/
const int NumberOfSteps = 8e2;
const double dstep=1./double(NumberOfSteps);


double bessel_polynomial(double rho, int l, int m){
    return std::cyl_bessel_j(double(m),double(bessel_zeros.at(m).at(l-1)*rho/radius));
}

double normalization_factor(int l1, int m1, int l2, int m2){
    double norm_fac = 1/(
        std::sqrt(2*pi*length*radius*radius*pow(std::cyl_bessel_j(double(m1), bessel_zeros.at(m1+1).at(l1-1)),2)) *
        std::sqrt(2*pi*length*radius*radius*pow(std::cyl_bessel_j(double(m2), bessel_zeros.at(m2+1).at(l2-1)),2))
    );
    return norm_fac;
}

double function_H(int m1, int l1, int m2, int l2, double rho1, double rho2){
    return bessel_polynomial(rho1, l1, m1)*bessel_polynomial(rho2, l1, m1);
}

double function_A_1(int n1, int m1, int l1, int n2, int m2, int l2, double rho1, double rho2){
    return rho1*rho1*(
        rho2*rho2*(
            length*length*(bessel_zeros.at(m1).at(l1-1)*bessel_zeros.at(m1).at(l1-1) + bessel_zeros.at(m2).at(l2-1))*bessel_zeros.at(m2).at(l2-1)
            + 4*pi*pi*(n1*n1+n2*n2)*radius*radius
            )+
        length*length*(mass*mass-m1*m1)*radius*radius
    )+length*length*(mass*mass-m2*m2)*radius*radius*rho2*rho2;
}
double function_A_2(int n1, int m1, int l1, int n2, int m2, int l2, double rho1, double rho2){
    return rho1*rho1*(
        rho2*rho2*(
            length*length*(bessel_zeros.at(m1).at(l1-1)*bessel_zeros.at(m1).at(l1-1) + bessel_zeros.at(m2).at(l2-1))*bessel_zeros.at(m2).at(l2-1)
            + 4*pi*pi*(n1*n1+n2*n2)*radius*radius
            )+
        length*length*(mass*mass-m2*m2)*radius*radius
    )+length*length*(mass*mass-m1*m1)*radius*radius*rho2*rho2;
}



double combined_integral(double g, int n1, int m1, int l1, int n2, int m2, int l2){
    double delta_integral = 0;
    double laplace_integral = 0;
    //double norm_fac = 1;
    double norm_fac = normalization_factor(l1, m1, l2, m2);
    if (n1==n2){
        for (int i0=0 ; i0<NumberOfSteps ; i0++){
            double rho1=i0*dstep*radius+1e-8;
            delta_integral += g*8*pi*length*dstep*radius*pow(function_H(m1, l1, m2, l2, rho1, rho1),2);
            for (int i1=0 ; i1<NumberOfSteps ; i1++){
                double rho2=i1*dstep*radius+1e-8;
                double firstH = function_H(m1, l1, m2, l2, rho2, rho1);
                double secondH = function_H(m1, l1, m2, l2, rho1, rho2);
                double firstA = function_A_1(n1, m1, l1, n2, m2, l2, rho1, rho2);
                double secondA = function_A_2(n1, m1, l1, n2, m2, l2, rho1, rho2);
                double laplace_integral_increment = -hbar*0.5/mass*4*pi*pi/radius/radius*dstep*radius*dstep*radius/rho1/rho1/rho2/rho2*(
                -firstH*firstH*firstA
                -secondH*secondH*secondA
                -firstH*secondH*firstA
                -firstH*secondH*secondA
                );
                laplace_integral += laplace_integral_increment;
            }
        }
    }
    else {
        for (int i0=0 ; i0<NumberOfSteps ; i0++){
            double rho1=i0*dstep*radius+1e-8;
            delta_integral += g*8*pi*length*dstep*radius*pow(function_H(m1, l1, m2, l2, rho1, rho1),2);
            for (int i1=0 ; i1<NumberOfSteps ; i1++){
                double rho2=i1*dstep*radius+1e-8;
                double firstH = function_H(m1, l1, m2, l2, rho2, rho1);
                double secondH = function_H(m1, l1, m2, l2, rho1, rho2);
                double firstA = function_A_1(n1, m1, l1, n2, m2, l2, rho1, rho2);
                double secondA = function_A_2(n1, m1, l1, n2, m2, l2, rho1, rho2);
                double laplace_integral_increment = -hbar*0.5/mass*4*pi*pi/radius/radius*dstep*radius*dstep*radius/rho1/rho1/rho2/rho2*(
                -firstH*firstH*firstA
                -secondH*secondH*secondA
                //-firstH*secondH*firstA
                //-firstH*secondH*secondA
                );
                laplace_integral += laplace_integral_increment;
            }
        }
    }
    return (laplace_integral+delta_integral)*norm_fac*norm_fac*0.5; //normalizacja bozonowa (1/2) i dla obu funkcji falowych (norm_fac^2)
}



int main(int argc, char** argv){
auto start = std::chrono::high_resolution_clock::now();

    double g=1;

    std::cout << "number of steps: " << NumberOfSteps << std::endl;
    //std::cout << bessel_zeros.at(0).at(1-1) << ' ' << bessel_zeros.at(5).at(3-1)<<std::endl;
    /*bessel_zeros.at(m).at(l-1)*/
/*
    std::cout <<
    8*pi*pi*(n1*n1+n2*n2)*laplace_integral(l1, l2, m)+2*pi*g*length*delta_integral(l1, l2, m)
    <<std::endl;
*/
    int nmin = 0;
    int nmax = 1;

    int lmin = 1;
    int lmax = 1;
    
    int mmin = 0;
    int mmax = 1;



    std::vector<std::string> parameters = {"n1", "n2", "l1", "l2", "m1", "m2", "integral_value"};
    std::vector<std::vector<double>> values;

    for (int n1 = nmin; n1 <= nmax; ++n1) {
        for (int n2 = nmin; n2 <= nmax; ++n2) {
            for (int l1 = lmin; l1 <= lmax; ++l1) {
                for (int l2 = lmin; l2 <= lmax; ++l2) {
                    for (int m1 = mmin; m1 <= mmax; ++m1) {
                        for (int m2 = mmin; m2 <= mmax; ++m2) {
                            std::cout << "n1: " << n1 << ", n2: " << n2 << ", l1: " << l1 
                                    << ", l2: " << l2 << ", m1: " << m1 << ", m2: "<< m2 << ", element tensora: ";

                            
                            double tensor_element = combined_integral(g, n1, m1, l1, n2, m2, l2);
                            /*
                            double tensor_element = 8*pi*pi*(n1*n1+n2*n2)*laplace_integral(l1, l2, m)+2*pi*g*length*delta_integral(l1, l2, m);
                            */
                            std::cout << tensor_element << std::endl;

                            values.push_back(std::vector<double>({double(n1),double(n2),double(l1),double(l2),double(m1),double(m2), tensor_element}));
                        }
                    }
                }
            }
        }
    }

auto end = std::chrono::high_resolution_clock::now();
double elapsed_seconds = std::chrono::duration<double>(end - start).count();
std::cout << "Elapsed time: " << elapsed_seconds << " seconds,\t g = " << g << std::endl;


std::ostringstream oss;
oss << "Elapsed time: " << elapsed_seconds << "\tg = " << g;
std::string extrainfo = oss.str();
    functions Function;
    Function.write_to_file("output.txt", parameters, values, extrainfo);

    return 0;
}
