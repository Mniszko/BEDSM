#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include "eigen/Eigen/Dense"
#include "./src/variables.h"
#include <chrono>
#include <fstream>
#include <string>
#include <sstream>

const long double hbar = 1; /*1.0545718e-34;*/
const double pi = 3.14159265358979323846;
const double sqrtPi = 1.7724538509055159;

const double length = 1; /*długość małego walca == L*/
const double radius = 1; /*promień małego walca == sigma*/
const double mass = 1;

const int NumberOfSteps = 8e2;
const double dstep=1./double(NumberOfSteps);


//równanie pod całkami radialnymi jest postaci:
/*
-(1/(\[Pi]^2 \[Sigma]^2))4 (L^2 ((Subscript[j, Subscript[m, 1],Subscript[l, 1]])^2+(Subscript[j, Subscript[m, 2],Subscript[l, 2]])^2)+4 \[Pi]^2 (\!\(
\*SubsuperscriptBox[\(n\), \(1\), \(2\)] + 
\*SubsuperscriptBox[\(n\), \(2\), \(2\)]\)) \[Sigma]^2) ((2 sin^2(\[Pi] (Subscript[m, 1]-Subscript[m, 2])) sin^2(\[Pi] (Subscript[n, 1]-Subscript[n, 2])) Subscript[J, Subscript[m, 1]]((Subscript[j, Subscript[m, 1],Subscript[l, 1]] Subscript[\[Rho], 1])/\[Sigma]) Subscript[J, Subscript[m, 1]]((Subscript[j, Subscript[m, 1],Subscript[l, 1]] Subscript[\[Rho], 2])/\[Sigma]) Subscript[J, Subscript[m, 2]]((Subscript[j, Subscript[m, 2],Subscript[l, 2]] Subscript[\[Rho], 1])/\[Sigma]) Subscript[J, Subscript[m, 2]]((Subscript[j, Subscript[m, 2],Subscript[l, 2]] Subscript[\[Rho], 2])/\[Sigma]))/((Subscript[m, 1]-Subscript[m, 2])^2 (Subscript[n, 1]-Subscript[n, 2])^2)+\[Pi]^4 (Subscript[J, Subscript[m, 1]]((Subscript[j, Subscript[m, 1],Subscript[l, 1]] Subscript[\[Rho], 2])/\[Sigma])^2 Subscript[J, Subscript[m, 2]]((Subscript[j, Subscript[m, 2],Subscript[l, 2]] Subscript[\[Rho], 1])/\[Sigma])^2+Subscript[J, Subscript[m, 1]]((Subscript[j, Subscript[m, 1],Subscript[l, 1]] Subscript[\[Rho], 1])/\[Sigma])^2 Subscript[J, Subscript[m, 2]]((Subscript[j, Subscript[m, 2],Subscript[l, 2]] Subscript[\[Rho], 2])/\[Sigma])^2))
*/


double bessel_polynomial(double rho, int l, int m){
    return std::cyl_bessel_j(double(m),double(bessel_zeros.at(m).at(l-1)*rho/radius));
}

double normalization_factor(int l1, int m1, int l2, int m2){
    double norm_fac = 1/(
        std::sqrt(2*pi*length*radius*radius*pow(std::cyl_bessel_j(double(m1+1), bessel_zeros.at(m1).at(l1-1)),2)) *
        std::sqrt(2*pi*length*radius*radius*pow(std::cyl_bessel_j(double(m2+1), bessel_zeros.at(m2).at(l2-1)),2))
    );
    return norm_fac;
}

double function_H(int m1, int l1, int m2, int l2, double rho1, double rho2){
    // ta funkcja jest niepoprawnie zdefiniowana (albo niepoprawnie wołana). Może problem jest z funkcją bessel_polynomial samą w sobie.
    return bessel_polynomial(rho1, l1, m1)*bessel_polynomial(rho2, l2, m2);
}

double function_A(int n1, int m1, int l1, int n2, int m2, int l2, double rho1, double rho2){
    return length*length*(bessel_zeros.at(m1).at(l1-1)*bessel_zeros.at(m1).at(l1-1) + bessel_zeros.at(m2).at(l2-1)*bessel_zeros.at(m2).at(l2-1))+4*pi*pi*(n1*n1+n2*n2)*radius*radius;
}

//tylko jeśli n1 != n2 i m1 != m2
double function_E(int n1, int m1, int l1, int n2, int m2, int l2, double rho1, double rho2){
    if (n1 != n2 && m1 != m2){
        return 2*std::sin(pi*(n1-n2))*std::sin(pi*(n1-n2))*
        std::sin(pi*(m1-m2))*std::sin(pi*(m1-m2))/((m1-m2)*(m1-m2)*(n1-n2)*(n1-n2));
    }
    else if (n1 != n2) {
        return 2*std::sin(pi*(n1-n2))*std::sin(pi*(n1-n2))/((n1-n2)*(n1-n2));
    }
    else if (m1 != m2) {
        return 2*std::sin(pi*(m1-m2))*std::sin(pi*(m1-m2))/((m1-m2)*(m1-m2));
    }
    else {
        return 2;
    }
}

double combined_integral(double g, int n1, int m1, int l1, int n2, int m2, int l2){
    double delta_integral = 0;
    double laplace_integral = 0;
    //double norm_fac = 1;
    double norm_fac = normalization_factor(l1, m1, l2, m2);
    if (n1==n2){
        for (int i0=0 ; i0<NumberOfSteps ; i0++){
            double rho1=i0*dstep*radius+1e-8;
            delta_integral += g*8*pi*length*dstep*radius*pow(function_H(m1, l1, m2, l2, rho1, rho1),2); //należy to policzyć ponownie
            for (int i1=0 ; i1<NumberOfSteps ; i1++){
                double rho2=i1*dstep*radius+1e-8;
                double funH2 = function_H(m1, l1, m2, l2, rho2, rho1);
                double funH1 = function_H(m1, l1, m2, l2, rho1, rho2);
                double funA = function_A(n1, m1, l1, n2, m2, l2, rho1, rho2);
                double funE = function_E(n1, m1, l1, n2, m2, l2, rho1, rho2);
                double laplace_integral_increment = hbar*0.5/mass/(pi*pi)/(radius*radius)*4*funA*dstep*radius*dstep*radius*(
                    funE*funH1*funH2+
                    pi*pi*pi*pi*(funH1*funH1+funH2*funH2));
                laplace_integral += laplace_integral_increment;
            }
        }
    }
    return (laplace_integral+delta_integral)*norm_fac*norm_fac*0.5; //normalizacja bozonowa (1/2) i dla obu funkcji falowych (norm_fac^2)
}

int main(int argc, char** argv){
auto start = std::chrono::high_resolution_clock::now();

    double g=0;

    std::cout << "number of steps: " << NumberOfSteps << std::endl;

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

//problemy są bezpośrednio w funkcji Bessela, spróbuję wyplotować ją tutaj i w mathematice żeby sprawdzić czy zadziała
// jest możliwość, że korzystam ze złych funkcji bessela (wina biblioteki?).