/*
program trzeba kompilować w
c++ 17 albo wyższym np. za pomocą
komendy:
"g++ -std=c++17 LL-GP_bessel.cpp"

----------------------------------

zastosowałem metodę małych kroków,
do kwantowych symulacji lepsza 
jest podobno Numerova, zwłaszcza,
że pozwala wierniej odwzorować
deltę (tutaj w. w. prostokąt o 
skończonej objętości)

----------------------------------

macierz energii jest symetryczna,
więc liczenie elementów 
pozadiagonalnych można łatwo 
przyspieszyć 2-krotnie

*/

#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include "eigen/Eigen/Dense" /*eigen znajduje się w folderze z programem*/


const double hbar = 1; /*1.0545718e-34;*/
const double pi = 3.14159265358979323846;
const double sqrtPi = 1.7724538509055159;

const double length = 1; /*długość małego walca == L*/
const double radius = 1; /*promień małego walca == sigma*/
const double mass = 1;

const int NumberOfSteps = 1e2;

/*rzędy dla kolejnych besseli, kolumny dla kolejnych zer danego bessela, tablica 10x10 przeklejona z obliczeń w mathematice*/
std::vector<std::vector<double>> bessel_zeros = {{2.40482555769577, 5.52007811028631, 8.65372791291101, 
  11.7915344390143, 14.9309177084878, 18.0710639679109, 
  21.2116366298793, 24.3524715307493, 27.4934791320403, 
  30.6346064684320}, {3.83170597020751, 7.01558666981560, 
  10.1734681350627, 13.3236919363142, 16.4706300508776, 
  19.6158585104682, 22.7600843805928, 25.9036720876184, 
  29.0468285349169, 32.1896799109744}, {5.13562230184068, 
  8.41724414039986, 11.6198411721491, 14.7959517823513, 
  17.9598194949878, 21.1169970530218, 24.2701123135731, 
  27.4205735499846, 30.5692044955164, 
  33.7165195092227}, {6.38016189592398, 9.76102312998172, 
  13.0152007216984, 16.2234661603188, 19.4094152264350, 
  22.5827295931044, 25.7481666992950, 28.9083507809218, 
  32.0648524070977, 35.2186707386101}, {7.58834243450398, 
  11.0647094885010, 14.3725366716176, 17.6159660498048, 
  20.8269329569624, 24.0190195247711, 27.1990877659813, 
  30.3710076671172, 33.5371377118192, 
  36.6990011287446}, {8.77148381596118, 12.3386041974669, 
  15.7001740797117, 18.9801338751799, 22.2177998965613, 
  25.4303411542227, 28.6266183072911, 31.8117167240478, 
  34.9887812945593, 38.1598685619671}, {9.93610952422099, 
  13.5892901705412, 17.0038196678160, 20.3207892135665, 
  23.5860844355814, 26.8201519834114, 30.0337223865705, 
  33.2330417628471, 36.4220196682585, 
  39.6032394160754}, {11.0863700192496, 14.8212687270132, 
  18.2875828324817, 21.6415410198484, 24.9349278876730, 
  28.1911884594832, 31.4227941922656, 34.6370893520693, 
  37.8387173828536, 41.0307736915855}, {12.2250922639885, 
  16.0377741908877, 19.5545364309971, 22.9451731318746, 
  26.2668146411766, 29.5456596709985, 32.7958000373415, 
  36.0256150638696, 39.2404479951781, 
  42.4438877432736}, {13.3543004774407, 17.2412203824891, 
  20.8070477892641, 24.2338852577506, 27.5837489635730, 
  30.8853789676967, 34.1543779238551, 37.4000999771566, 
  40.6285537189645, 43.8438014203373}, {14.4755006865580, 
  18.4334636669666, 22.0469853646978, 25.5094505541828, 
  28.8873750635305, 32.2118561997127, 35.4999092053739, 
  38.7618070178817, 42.0041902366718, 45.2315741035350}};

std::complex<double> Z_component(double z, int n){
    return std::exp(std::complex<double>(0,2*z*n*pi/length))/std::sqrt(length);
}
std::complex<double> Phi_component(double phi) {
    return std::exp(std::complex<double>(0,mass*phi));
}
std::complex<double> Rho_component(double rho, int l, int m){
    //std::cout << double(bessel_zeros.at(l).at(m)*rho/radius) <<std::endl;
    //std::cout << rho <<std::endl;
    return std::cyl_bessel_j(double(l),double(bessel_zeros.at(l).at(m)*rho/radius));
}

std::complex<double> Psi_symmetrized(double rho1, double phi1, double z1, int n1, int l1, int m1, double rho2, double phi2, double z2, int n2, int l2, int m2){
    return Z_component(z1, n1)*Phi_component(phi1)*Rho_component(rho1, l1, m1)*
        Z_component(z2, n2)*Phi_component(phi2)*Rho_component(rho2, l2, m2) +
        Z_component(z2, n1)*Phi_component(phi2)*Rho_component(rho2, l1, m1)*
        Z_component(z1, n2)*Phi_component(phi1)*Rho_component(rho1, l2, m2);
}

/*przyjąłem zakres działania delty równy dwukrotności kroku*/
std::complex<double> delta(double rho1, double phi1, double z1, double rho2, double phi2, double z2){
    if (abs(rho1-rho2)<=2./NumberOfSteps*length && abs(phi1-phi2)<=4./NumberOfSteps*pi && abs(z1-z2)<=2./NumberOfSteps*length){
        return std::complex<double>(0.5*NumberOfSteps/sqrtPi,0);
    }
    else {return 0;}
}

std::complex<double> hamiltonian_times_psi(double rho1, double phi1, double z1, int n1, int l1, int m1, double rho2, double phi2, double z2, int n2, int l2, int m2, double g){
    double h = 1e-9; /*niezależne od kroku*/
    std::complex<double> DoublePsi = std::complex<double> (2,0)*Psi_symmetrized(rho1, phi1, z1, n1, l1, m1, rho2, phi2, z2, n2, l2, m2);
    return
        hbar*hbar/2/mass*( /*prawdopodobnie można to przyspieszyć ok. 2 krotnie pomijając zmiany rho2, bo jest efektywnie identyczne z rho1, być może da się pominąć complex double*/
        /*pierwsze pochodne po promieniu*/
            (Psi_symmetrized(rho1+h*radius, phi1, z1, n1, l1, m1, rho2, phi2, z2, n2, l2, m2)-Psi_symmetrized(rho1-h*radius, phi1, z1, n1, l1, m1, rho2, phi2, z2, n2, l2, m2))*std::complex<double>(1/h/radius/rho1,0) + 
            (Psi_symmetrized(rho1+h*radius, phi1, z1, n1, l1, m1, rho2, phi2, z2, n2, l2, m2)-Psi_symmetrized(rho1-h*radius, phi1, z1, n1, l1, m1, rho2, phi2, z2, n2, l2, m2))*std::complex<double>(1/h/radius/rho2,0) + 
        /*druga pochodna po promieniu*/
            (Psi_symmetrized(rho1+h*radius, phi1, z1, n1, l1, m1, rho2, phi2, z2, n2, l2, m2)-
DoublePsi+Psi_symmetrized(rho1-h*radius, phi1, z1, n1, l1, m1, rho2, phi2, z2, n2, l2, m2))*std::complex<double>(1/h/h/radius/radius,0) + 
            (Psi_symmetrized(rho1, phi1, z1, n1, l1, m1, rho2+h*radius, phi2, z2, n2, l2, m2)-
DoublePsi+Psi_symmetrized(rho1, phi1, z1, n1, l1, m1, rho2-h*radius, phi2, z2, n2, l2, m2))*std::complex<double>(1/h/h/radius/radius,0) +
        /*druga pochodna po kącie*/
            (Psi_symmetrized(rho1, phi1+h*2*pi, z1, n1, l1, m1, rho2, phi2, z2, n2, l2, m2)-
DoublePsi+Psi_symmetrized(rho1, phi1, z1, n1, l1, m1, rho2, phi2-h*2*pi, z2, n2, l2, m2))*std::complex<double>(1/h/h/pi/pi/rho1/rho1,0) + 
            (Psi_symmetrized(rho1, phi1, z1, n1, l1, m1, rho2, phi2+h*2*pi, z2, n2, l2, m2)-
DoublePsi+Psi_symmetrized(rho1, phi1, z1, n1, l1, m1, rho2, phi2-h*2*pi, z2, n2, l2, m2))*std::complex<double>(1/h/h/pi/pi/rho2/rho2,0) +
        /*druga pochodna po z*/
            (Psi_symmetrized(rho1, phi1, z1+h*length, n1, l1, m1, rho2, phi2, z2, n2, l2, m2)-
DoublePsi+Psi_symmetrized(rho1, phi1, z1-h*length, n1, l1, m1, rho2, phi2, z2, n2, l2, m2))*std::complex<double>(1/h/h/length/length,0) +
            (Psi_symmetrized(rho1, phi1, z1, n1, l1, m1, rho2, phi2, z2+h*length, n2, l2, m2)-
DoublePsi+Psi_symmetrized(rho1, phi1, z1, n1, l1, m1, rho2, phi2, z2-h*length, n2, l2, m2))*std::complex<double>(1/h/h/length/length,0)
        ) + std::complex<double>(g/2,0)*DoublePsi*delta(rho1, phi1, z1, rho2, phi2, z2)
        ;
}

double single_integral (int n1, int l1, int m1, int n2, int l2, int m2, double g){
    double integrated_value = 0;
    for (int i0=0 ; i0<NumberOfSteps ; i0++){
            for (int i1=0 ; i1<NumberOfSteps ; i1++){
                for (int i2=0 ; i2<NumberOfSteps ; i2++){
                    for (int i3=0 ; i3<NumberOfSteps ; i3++){
                        for (int i4=0 ; i4<NumberOfSteps ; i4++){
                            for (int i5=0 ; i5<NumberOfSteps ; i5++){
                            double rho1=i0*radius+1e-8;
                            double phi1=i1*2*pi+1e-8;
                            double z1=i2*length+1e-8;
                            double rho2=i3*radius+1e-8;
                            double phi2=i4*2*pi+1e-8;
                            double z2=i5*length+1e-8;

                            integrated_value +=  abs(std::conj(
                                Psi_symmetrized(rho1, phi1, z1, n1, l1, m1, rho2, phi2, z2, n2, l2, m2)
                            )*hamiltonian_times_psi(rho1, phi1, z1, n1, l1, m1, rho2, phi2, z2, n2, l2, m2, g));
                            }     
                        }        
                    }        
                }        
            }       
    }
    return integrated_value/(radius*radius*4*pi*pi*length*length);
}



int main(int argc, char** argv){
    /*poniżej główna pętla numeryczna ze stałymi krokami całkowania, niezoptymalizowana*/
    std::cout << single_integral(1, 1, 1, 1, 1, 1, 2) << std::endl;
    return 0;

}
