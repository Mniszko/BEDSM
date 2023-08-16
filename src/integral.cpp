#include "integral.h"

#include <functional>
#include <iostream>
#include <stdexcept>
#include <gsl/gsl_integration.h>
#include <boost/math/special_functions/hypergeometric_1F1.hpp>


const long double hbar = 1; /*1.0545718e-34;*/
const double pi = 3.14159265358979323846;
const double sqrtPi = 1.7724538509055159;

//const double length = 1; /*długość małego walca == L*/
//const double rev_length = 1/length;
const double radius = 1; /*promień małego walca == sigma*/
const double mass = 1;

const int NumberOfSteps = 8e2;
const double dstep=1./double(NumberOfSteps);
const double h = 1;

// Simple RAII wrapper 
class IntegrationWorkspace {
    gsl_integration_workspace * wsp;

    public:
    IntegrationWorkspace(const size_t n=1000):
        wsp(gsl_integration_workspace_alloc(n)) {}
    ~IntegrationWorkspace() { gsl_integration_workspace_free(wsp); }

    operator gsl_integration_workspace*() { return wsp; }
};

// Build gsl_function from lambda
template <typename F>
class gsl_function_pp: public gsl_function {
    const F func;
    static double invoke(double x, void *params) {
        return static_cast<gsl_function_pp*>(params)->func(x);
    }
    public:
        gsl_function_pp(const F& f) : func(f) {
        function = &gsl_function_pp::invoke; //inherited from gsl_function
        params   = this;                     //inherited from gsl_function
        }
        operator gsl_function*(){return this;}
};

// Helper function for template construction
template <typename F>
gsl_function_pp<F> make_gsl_function(const F& func) {
    return gsl_function_pp<F>(func);
}

void CompleteIntegral::check_if_valid_key(int k){
    if(k != 1 && k != 2 && k != 3 && k != 4 && k != 5 && k != 6){
                throw std::runtime_error("Wrong value for gaussian interpolation key");
            }
}

void CompleteIntegral::print_result(){
            std::cout << "result:\t" << this->result << "\nerror:\t" << this->error << std::endl;
        }

void CompleteIntegral::changeMultiIndex(int n1, int m1, int u1, int n2, int m2, int u2){
    this->n1 = n1;
    this->m1 = m1;
    this->u1 = u1;
    this->n2 = n2;
    this->m2 = m2;
    this->u2 = u2;
}
void CompleteIntegral::changeLimit(int limit){
    this->limit = limit;
}

double CompleteIntegral::norm_psi(double r1, void* params) {
    double* params_arr = static_cast<double*>(params);
    double h1 = params_arr[0];
    double m1 = params_arr[1];
    double u1 = params_arr[2];

    double exponent = -h1 * r1 * r1;
    double r_power = pow(r1, 2 * m1);
    double hypergeometric1F1 = boost::math::hypergeometric_1F1(-u1, 1 + m1, h1 * r1 * r1);

    return exp(exponent) * r_power * hypergeometric1F1 * hypergeometric1F1;
}

void CompleteIntegral::integrate_norm() {
    double result1, error1, result2, error2;
    IntegrationWorkspace wsp1(this->limit);

    gsl_function F = make_gsl_function( [&](double rho) {
            // inside of a function
            double exponent = -h * rho * rho;
            double r_power = std::pow(rho, (double)2 * this->m1);
            double hypergeometric1F1 = boost::math::hypergeometric_1F1(-this->u1, 1 + this->m1, h * rho * rho);

            return exp(exponent) * r_power * hypergeometric1F1 * hypergeometric1F1; //jeśli tu pojawi się problem to całkę trzeba będzie zastąpić całką nie laplacianu - tylko oscylatora. 
        } );
    gsl_integration_qagiu(&F, 0, this->epsabs, this->epsrel, this->limit, wsp1, &result1, &error1);

    if (this->n1!=this->n2 || this->m1!=this->m2 || this->u1!=this->u2){
        IntegrationWorkspace wsp2(this->limit);
        gsl_function F = make_gsl_function( [&](double rho) {
            // inside of a function
            double exponent = -h * rho * rho;
            double r_power = pow(rho, 2 * this->m2);
            double hypergeometric1F1 = boost::math::hypergeometric_1F1(-this->u2, 1 + this->m2, h * rho * rho);

            return exp(exponent) * r_power * hypergeometric1F1 * hypergeometric1F1; //jeśli tu pojawi się problem to całkę trzeba będzie zastąpić całką nie laplacianu - tylko oscylatora. 
        } );
        gsl_integration_qagiu(&F, 0, this->epsabs, this->epsrel, this->limit, wsp2, &result2, &error2);
    } else {
        result2 = result1;
        error2 = error1;
    }
    this->norm1 = result1*2*pi;
    this->norm2 = result2*2*pi;
}

/*
double factorial_of(int x){
    int factorial;
    for(int i = 1; i <= x; ++i) {
            factorial *= i;
        }
    return (double)factorial;
}

void CompleteIntegral::integrate_norm(){
    double result1, result2;
    result1 = factorial_of(this->u1+this->m1)/factorial_of(this->u1)/factorial_of(this->m1)/factorial_of(this->m1);
    result2 = factorial_of(this->u2+this->m2)/factorial_of(this->u2)/factorial_of(this->m2)/factorial_of(this->m2);
    this->norm1 = result1*2*pi;
    this->norm2 = result2*2*pi;
}
*/

//definitions of different integrals
double CompleteIntegral::rhoPhiLaplacianIntegral_real(int m, int u){//integral of Psi*R
    double result, abserr, inner_result, inner_abserr;
    double * result_and_error = new double[2];
    IntegrationWorkspace wsp1(this->limit);
    IntegrationWorkspace wsp2(this->limit);
    auto outer = make_gsl_function( [&](double phi) {
        auto inner = make_gsl_function( [&](double rho) {
            // inside of a function
            double phi_component = std::cos(m*phi);
            double rho_component = std::exp(-h*rho*rho*0.5)*pow(rho,(double)m)*boost::math::hypergeometric_1F1(-u, m + 1, rho*rho * h);
            double real_psi = phi_component*rho_component;
            double real_laplacianPsi = (2*u+m+1-rho*rho)*phi_component*rho_component;
            return real_psi*real_laplacianPsi*this->length; //jeśli tu pojawi się problem to całkę trzeba będzie zastąpić całką nie laplacianu - tylko oscylatora. 
        } ); //mapped to number of gaussian quadrature points
        gsl_integration_qagiu(inner, 0, epsabs, epsrel, limit, wsp1,
                                &inner_result, &inner_abserr);
        return inner_result;
        });
    int key=this->key;
    gsl_integration_qag(outer, 0, 2*pi, epsabs, epsrel, limit, key, wsp1,&result, &abserr);
    result_and_error[0] = result;
    result_and_error[1] = abserr;
    this->result += result;
    this->error += abserr;
    return result;
}
double CompleteIntegral::rhoPhiLaplacianIntegral_imaginary(int m, int u){//integral of Psi*R
    double result, abserr, inner_result, inner_abserr;
    double * result_and_error = new double[2];
    IntegrationWorkspace wsp1(this->limit);
    IntegrationWorkspace wsp2(this->limit);
    auto outer = make_gsl_function( [&](double phi) {
        auto inner = make_gsl_function( [&](double rho) {
            // inside of a function
            double phi_component = std::sin(m*phi);
            double rho_component = std::exp(-h*rho*rho*0.5)*pow(rho,(double)m)*boost::math::hypergeometric_1F1(-u, m + 1, rho*rho * h);
            double imaginary_psi = phi_component*rho_component;
            double imaginary_laplacianPsi = (2*u+m+1-rho*rho)*imaginary_psi;
            return imaginary_psi*imaginary_laplacianPsi*this->length; //jeśli tu pojawi się problem to całkę trzeba będzie zastąpić całką nie laplacianu - tylko oscylatora. 
        } ); //mapped to number of gaussian quadrature points
        gsl_integration_qagiu(inner, 0, epsabs, epsrel, limit, wsp1,
                                &inner_result, &inner_abserr);
        return inner_result;
        });
    int key=this->key;
    gsl_integration_qag(outer, 0, 2*pi, epsabs, epsrel, limit, key, wsp1,&result, &abserr);
    result_and_error[0] = result;
    result_and_error[1] = abserr;
    this->result += result;
    this->error += abserr;
    return result;
}
double CompleteIntegral::rhoPhiLaplacianIntegralWithHarmonic_real(int m, int u){//integral of Psi*R
    double result, abserr, inner_result, inner_abserr;
    double * result_and_error = new double[2];
    IntegrationWorkspace wsp1(this->limit);
    IntegrationWorkspace wsp2(this->limit);
    auto outer = make_gsl_function( [&](double phi) {
        auto inner = make_gsl_function( [&](double rho) {
            // inside of a function
            double phi_component = std::cos(m*phi);
            double rho_component = std::exp(-h*rho*rho*0.5)*pow(rho,(double)m)*boost::math::hypergeometric_1F1(-u, m + 1, rho*rho * h);
            double real_psi = phi_component*rho_component;
            double real_laplacianPsi = (2*u+m+1)*phi_component*rho_component;
            return real_psi*real_laplacianPsi*this->length; //jeśli tu pojawi się problem to całkę trzeba będzie zastąpić całką nie laplacianu - tylko oscylatora. 
        } ); //mapped to number of gaussian quadrature points
        gsl_integration_qagiu(inner, 0, epsabs, epsrel, limit, wsp1,
                                &inner_result, &inner_abserr);
        return inner_result;
        });
    int key=this->key;
    gsl_integration_qag(outer, 0, 2*pi, epsabs, epsrel, limit, key, wsp1,&result, &abserr);
    result_and_error[0] = result;
    result_and_error[1] = abserr;
    this->result += result;
    this->error += abserr;
    return result;
}
double CompleteIntegral::rhoPhiLaplacianIntegralWithHarmonic_imaginary(int m, int u){//integral of Psi*R
    double result, abserr, inner_result, inner_abserr;
    double * result_and_error = new double[2];
    IntegrationWorkspace wsp1(this->limit);
    IntegrationWorkspace wsp2(this->limit);
    auto outer = make_gsl_function( [&](double phi) {
        auto inner = make_gsl_function( [&](double rho) {
            // inside of a function
            double phi_component = std::sin(m*phi);
            double rho_component = std::exp(-h*rho*rho*0.5)*pow(rho,(double)m)*boost::math::hypergeometric_1F1(-u, m + 1, rho*rho * h);
            double imaginary_psi = phi_component*rho_component;
            double imaginary_laplacianPsi = (2*u+m+1)*imaginary_psi;
            return imaginary_psi*imaginary_laplacianPsi*this->length; //jeśli tu pojawi się problem to całkę trzeba będzie zastąpić całką nie laplacianu - tylko oscylatora. 
        } ); //mapped to number of gaussian quadrature points
        gsl_integration_qagiu(inner, 0, epsabs, epsrel, limit, wsp1,
                                &inner_result, &inner_abserr);
        return inner_result;
        });
    int key=this->key;
    gsl_integration_qag(outer, 0, 2*pi, epsabs, epsrel, limit, key, wsp1,&result, &abserr);
    result_and_error[0] = result;
    result_and_error[1] = abserr;
    this->result += result;
    this->error += abserr;
    return result;
}

//you could just rewrite that as single integral * 2 pi, it wouldn't change anything imo
double CompleteIntegral::rhoPhiHarmonicPotentialIntegral(int m, int u){//integral of Psi*R
    double result, abserr, inner_result, inner_abserr;
    double * result_and_error = new double[2];
    IntegrationWorkspace wsp1(this->limit);
    IntegrationWorkspace wsp2(this->limit);
    auto outer = make_gsl_function( [&](double phi) {
        auto inner = make_gsl_function( [&](double rho) {
            // inside of a function
            double psi = std::exp(-h*rho*rho*0.5)*pow(rho,(double)m)*boost::math::hypergeometric_1F1(-u, m + 1, rho*rho * h);
            double potentialPsi = (rho*rho)*psi;
            return psi*potentialPsi*this->length; //this shouldn't be a problem because kinetic energy should remain in numerical boundaries (and there are no singularities)
        } ); //mapped to number of gaussian quadrature points
        gsl_integration_qagiu(inner, 0, epsabs, epsrel, limit, wsp1,
                                &inner_result, &inner_abserr);
        return inner_result;
        });
    int key=this->key;
    gsl_integration_qag(outer, 0, 2*pi, epsabs, epsrel, limit, key, wsp1,&result, &abserr);
    result_and_error[0] = result;
    result_and_error[1] = abserr;
    this->result += result;
    this->error += abserr;
    return result;
}

double CompleteIntegral::zLaplacianIntegral_real(int n, double norm){//integral of Z
    double result, abserr;
    double * result_and_error = new double[2];
    IntegrationWorkspace wsp(this->limit);
    auto integrand = make_gsl_function( [&](double z) {
        // inside of a function
            double real_z = cos(2*pi*n*z*this->rev_length);
            double real_laplacianZ = 4*pi*pi*this->rev_length*this->rev_length*n*n*real_z;
            return real_z*real_laplacianZ*norm;
        } );
    int key=this->key;
    gsl_integration_qag(integrand, 0, this->length, epsabs, epsrel, limit, key, wsp,&result, &abserr);
    result_and_error[0] = result;
    result_and_error[1] = abserr;
    this->result += result;
    this->error += abserr;
    return result;
}
double CompleteIntegral::zLaplacianIntegral_imaginary(int n, double norm){//integral of Z
    double result, abserr;
    double * result_and_error = new double[2];
    IntegrationWorkspace wsp(this->limit);
    auto integrand = make_gsl_function( [&](double z) {
        // inside of a function
            double imaginary_z = sin(2*pi*n*z*this->rev_length);
            double imaginary_laplacianZ = 4*pi*pi*this->rev_length*this->rev_length*n*n*imaginary_z;
            return imaginary_z*imaginary_laplacianZ*norm;
        } );
    int key=this->key;
    gsl_integration_qag(integrand, 0, this->length, epsabs, epsrel, limit, key, wsp,&result, &abserr);
    result_and_error[0] = result;
    result_and_error[1] = abserr;
    this->result += result;
    this->error += abserr;
    return result;
}

double CompleteIntegral::fastLaplacianWithHarmonicEnergy(){
    double result = (this->n1*this->n1+this->n2*this->n2)*4*pi*pi*this->rev_length*this->rev_length + 2*this->u1+2*this->u2+this->m1+this->m2+2;
    this->result += result;
    return result;
}

double CompleteIntegral::fast_add_over_harmonic(){
    double result;
    result += this->fastLaplacianWithHarmonicEnergy();
    return result;
}

double CompleteIntegral::test_integral(){
    double result = 0;
    result += this->zLaplacianIntegral_real(this->n2, this->norm2);
    result += this->zLaplacianIntegral_imaginary(this->n2, this->norm2);
    result += this->rhoPhiLaplacianIntegral_real(this->m2, this->u2);
    result += this->rhoPhiLaplacianIntegral_imaginary(this->m2, this->u2);
    result += this->rhoPhiHarmonicPotentialIntegral(this->m2,this->u2);
    result /= this->norm2/this->length; // Nie mogę spodziewać się integerów, trzeba pamiętać, że zIntegral daje energie będące wielokrotnościami PI. Żeby to  naprawić mógłbym przyjąć długość odpowiednio zależną od PI
    return result;
}

double CompleteIntegral::Integrate_over_harmonic(){ //troubleshooting and testing to be done
    double same_index_result=0, one_index_result=0, different_index_result=0;
    if(n1==n2 && m1==m2 && u1==u2){
        same_index_result += this->zLaplacianIntegral_real(this->n1, this->norm1);
        same_index_result += this->zLaplacianIntegral_imaginary(this->n1, this->norm1);
        same_index_result += this->rhoPhiLaplacianIntegral_real(this->m1, this->u1);
        same_index_result += this->rhoPhiLaplacianIntegral_imaginary(this->m1, this->u1);
        same_index_result += this->rhoPhiHarmonicPotentialIntegral(this->m1,this->u1);
        same_index_result *= 4./this->norm1*this->rev_length*0.5; //dzielenie przez dwa w związku z normalizacją, mnożenie przez osiem w związku z rozbiciem całki po symetryzacji
    }
    else{
        one_index_result += this->zLaplacianIntegral_real(this->n1, this->norm1);
        one_index_result += this->zLaplacianIntegral_imaginary(this->n1, this->norm1);
        one_index_result += this->rhoPhiLaplacianIntegral_real(this->m1, this->u1);
        one_index_result += this->rhoPhiLaplacianIntegral_imaginary(this->m1, this->u1);
        one_index_result += this->rhoPhiHarmonicPotentialIntegral(this->m1,this->u1);
        one_index_result *= 2./this->norm1*this->rev_length*0.5;

        different_index_result += this->zLaplacianIntegral_real(this->n2, this->norm2);
        different_index_result += this->zLaplacianIntegral_imaginary(this->n2, this->norm2);
        different_index_result += this->rhoPhiLaplacianIntegral_real(this->m2, this->u2);
        different_index_result += this->rhoPhiLaplacianIntegral_imaginary(this->m2, this->u2);
        different_index_result += this->rhoPhiHarmonicPotentialIntegral(this->m2,this->u2);
        different_index_result *= 2./this->norm2*this->rev_length*0.5;
    }
    double* return_value = new double[2];
    return_value[0] = this->result; 
    return_value[1] = this->error;
    return same_index_result+one_index_result+different_index_result;
}

double CompleteIntegral::Integrate_over_harmonic_with_delta_potential_1D_motion_head_crush(double g){
    double result;
    if (abs(this->n1) != abs(this->n2)){
        result = this->Integrate_over_delta_1D_movement_real(g);
    } else {
        result = this->Integrate_over_delta_1D_movement_real(g);
        result += this->Integrate_over_harmonic();
    }
    return result;
}

double CompleteIntegral::Integrate_over_harmonic_with_delta(){
    return (double)0;
}//3D movement case without magnetic moment interactions
double CompleteIntegral::Integrate_over_delta_1D_movement_real(double g){
    double norm = sqrt(1/this->norm1/this->norm2);
    double result, abserr, middle_result, middle_abserr, inner_result, inner_abserr;
    IntegrationWorkspace wsp1(this->limit);
    IntegrationWorkspace wsp2(this->limit);
    IntegrationWorkspace wsp3(this->limit);
    auto outer = make_gsl_function([&](double z){
        auto middle = make_gsl_function( [&](double phi) {
            auto inner = make_gsl_function( [&](double rho) {
                // inside of a function
                double zpart1 = sin(2*pi*this->n1*z*this->rev_length);
                double phipart1 = sin(this->m1*phi);
                double rhopart1 = std::exp(-h*rho*rho*0.5)*pow(rho,(double)this->m1)*boost::math::hypergeometric_1F1(-this->u1, this->m1 + 1, rho*rho * h);
                double zpart2 = sin(2*pi*this->n2*z*this->rev_length);
                double phipart2 = sin(this->m2*phi);
                double rhopart2 = std::exp(-h*rho*rho*0.5)*pow(rho,(double)this->m2)*boost::math::hypergeometric_1F1(-this->u2, this->m2 + 1, rho*rho * h);
                double psi1 = rhopart1;
                double psi2 = rhopart2;
                return g*psi1*psi1*psi2*psi2*norm*norm*2*this->rev_length*this->rev_length; 
            } ); //mapped to number of gaussian quadrature points
            gsl_integration_qagiu(inner, 0, epsabs, epsrel, limit, wsp1,
                                    &inner_result, &inner_abserr);
            return inner_result;
            });
        int key=this->key;
        gsl_integration_qag(middle, 0, 2*pi, epsabs, epsrel, limit, key, wsp2,&middle_result, &middle_abserr);
        return middle_result;
    });
    int key=this->key;
    gsl_integration_qag(outer, 0, this->length, epsabs, epsrel, limit, key, wsp3, &result, &abserr);
    this->result += result;
    this->error += abserr;
    return result;
}

double CompleteIntegral::Integrate_over_delta_1D_movement_imaginary(double g){
    double norm = sqrt(1/this->norm1/this->norm2);
    double result, abserr, middle_result, middle_abserr, inner_result, inner_abserr;
    IntegrationWorkspace wsp1(this->limit);
    IntegrationWorkspace wsp2(this->limit);
    IntegrationWorkspace wsp3(this->limit);
    auto outer = make_gsl_function([&](double z){
        auto middle = make_gsl_function( [&](double phi) {
            auto inner = make_gsl_function( [&](double rho) {
                // inside of a function
                double zpart1 = sin(2*pi*this->n1*z*this->rev_length);
                double phipart1 = sin(this->m1*phi);
                double rhopart1 = std::exp(-h*rho*rho*0.5)*pow(rho,(double)this->m1)*boost::math::hypergeometric_1F1(-this->u1, this->m1 + 1, rho*rho * h);
                double zpart2 = sin(2*pi*this->n2*z*this->rev_length);
                double phipart2 = sin(this->m2*phi);
                double rhopart2 = std::exp(-h*rho*rho*0.5)*pow(rho,(double)this->m2)*boost::math::hypergeometric_1F1(-this->u2, this->m2 + 1, rho*rho * h);
                double psi1 = zpart1*phipart1*rhopart1;
                double psi2 = zpart2*phipart2*rhopart2;
                return g*psi1*psi1*psi2*psi2*norm; 
            } ); //mapped to number of gaussian quadrature points
            gsl_integration_qagiu(inner, 0, epsabs, epsrel, limit, wsp1,
                                    &inner_result, &inner_abserr);
            return inner_result;
            });
        int key=this->key;
        gsl_integration_qag(middle, 0, 2*pi, epsabs, epsrel, limit, key, wsp2,&middle_result, &middle_abserr);
        return middle_result;
    });
    int key=this->key;
    gsl_integration_qag(outer, 0, this->length, epsabs, epsrel, limit, key, wsp3, &result, &abserr);
    this->result += result;
    this->error += abserr;
    return result;
}



double* CompleteIntegral::Integrate_over_harmonic_with_delta_1D_movement(double g){//returns both diagonal elements of an array and delta integrated array (symmetric one)
    if(m1!=m2 || u1!=u2){
        std::cout << "Warning! - there should be only one value m and u for this function!" << std::endl;
    }
    double harmonic_result=0, delta_result=0;
    double* combined_result = new double[2];

    harmonic_result = this->Integrate_over_harmonic();
    delta_result = this->Integrate_over_delta_1D_movement_real(g);
    delta_result += this->Integrate_over_delta_1D_movement_imaginary(g);

    combined_result[0]=harmonic_result;
    combined_result[1]=delta_result;
    return combined_result;
}//replication of Lieb-Lienieger energy result in relation to variable g
