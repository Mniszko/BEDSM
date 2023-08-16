// following code is intended to run after call from python file "matrix.py", although it can be used as a standalone

#include "./src/integral.h" // stores integral class with all functionalities
#include "./src/variables.h" // besselj zeros aren't currently needed, stores writing to file function
#include "./src/flag_parser.h" // used to parse flags
#include <iostream>
#include <cmath>
#include <chrono>
#include <sstream>
#include <exception>
#include <string>
#include <gsl/gsl_errno.h>

//handling gnu type exceptions (handy to define here)
class GSLException : public std::exception {
private:
    std::string message;
public:
    GSLException(const std::string& msg) : message(msg) {}

    virtual const char* what() const throw() {
        return message.c_str();
    }
};


void custom_gsl_error_handler(const char * reason, const char * file, int line, int gsl_errno) {
    throw GSLException(std::string(reason) + " in file " + file + " at line " + std::to_string(line));
}



int main(int argc, char* argv[]){
    double g = 1;
    double length = 1;
    double epsabs = 1e-8;
    double epsrel = 1e-6;
    size_t limit = 100;
    int nmin = 0;
    int nmax = 2;

    int umin = 0;
    int umax = 0;
    
    int mmin = 0;
    int mmax = 0;
    
    FlagParser flag_parser(argc, argv);
    flag_parser.parse_flags();

    if (flag_parser.change_lFlag){
        length = flag_parser.change_lValue;
    }
    if (flag_parser.change_gFlag){
        g = flag_parser.change_gValue;
    }    
    if (flag_parser.minmaxFlag){
            nmin = flag_parser.nmin;
            nmax = flag_parser.nmax;
            mmin = flag_parser.mmin;
            mmax = flag_parser.mmax;
            umin = flag_parser.umin;
            umax = flag_parser.umax;
    }
    
    gsl_set_error_handler(custom_gsl_error_handler);

    auto start = std::chrono::high_resolution_clock::now();
    
    
    int n1 = 0;
    int m1 = 0;
    int u1 = 0;
    int n2 = 0;
    int m2 = 0;
    int u2 = 0;
    if (flag_parser.debugFlag && flag_parser.debugValues.size()==3){ //debugging purpose
        n1 = flag_parser.debugValues.at(0);
        m1 = flag_parser.debugValues.at(1);
        u1 = flag_parser.debugValues.at(2);
    } else if (flag_parser.debugFlag && flag_parser.debugValues.size()==6){ //debugging purpose
        n1 = flag_parser.debugValues.at(0);
        m1 = flag_parser.debugValues.at(1);
        u1 = flag_parser.debugValues.at(2);
        n2 = flag_parser.debugValues.at(3);
        m2 = flag_parser.debugValues.at(4);
        u2 = flag_parser.debugValues.at(5);
    } else {//main loop writing matrix in another file

            
            
            std::vector<std::string> parameters = {"n1", "n2", "u1", "u2", "m1", "m2", "integral_value"};
            std::vector<std::vector<double>> values;

            //initializing integral object
            CompleteIntegral integral(limit,epsabs,epsrel);
            integral.change_length(length);
            
            integral.change_key(6); //sets gaussian quadrature approximation resolution (integer 0<key<7) (strongly recommended to keep key value above 2)
            for (int n1 = nmin; n1 <= nmax; ++n1) {
            for (int n2 = nmin; n2 <= nmax; ++n2) { // opposite quantum numbers
                for (int u1 = umin; u1 <= umax; ++u1) {
                for (int u2 = umin; u2 <= umax; ++u2) {
                    for (int m1 = mmin; m1 <= mmax; ++m1) {
                    for (int m2 = mmin; m2 <= mmax; ++m2) {
                        //debugging purposes:
                        std::cout << "n1: " << n1 << ", n2: " << n2 << ", u1: " << u1 << ", u2: " << u2 << ", m1: " << m1 << ", m2: "<< m2 << ", element tensora: ";

                        //sets multiindex, integrates norm for both wavefunctions and stores it inside object
                        integral.changeMultiIndex(n1, m1, u1, n2, m2, u2);
                        integral.integrate_norm();
                        double tensor_element=0;

                        //try{
                        
                        tensor_element += integral.Integrate_over_harmonic();
                        //tensor_element += integral.fast_add_over_harmonic();
                        //tensor_element += integral.Integrate_over_delta_1D_movement_real(g);
                        
                        //tensor_element = integral.test_integral();
                        /*
                        } catch (const GSLException& e){
                            std::cout <<  e.what() << '\n' << "changing gaussian integration key" << std::endl;
                            m2 = m2-1;
                            integral.add_to_key();
                            continue;
                        }
                        */
                        //integral.change_key(1);
                        
                        //function responsible for integration.
                        //only delta potential
                        //integrate over harmonic saves wrong value as inner result of an integral (I may change that later), so one has to use external variable (here "tensor_element")

                        //integral.test_integral();
                        //integral.Integrate_over_delta_1D_movement_imaginary(g);
                        

                        //double tensor_element = integral.Integrate_over_harmonic_with_delta_potential_1D_motion_head_crush(g);
                        
                        //integral.Integrate_over_harmonic_with_delta_potential_1D_motion_head_crush(g);
                        //double tensor_element = integral.get_result();

                        //double tensor_element = integral.Integrate_over_harmonic_with_delta_1D_movement(g)[1];

                        //double tensor_element = integral.Integrate_over_harmonic();
                        
                        //resets error and result to 0
                        integral.restart_counter();
                        //integrates, applies normalization

                        //debugging purposes:
                        std::cout << tensor_element << std::endl;

                        values.push_back(std::vector<double>({double(n1),double(n2),double(u1),double(u2),double(m1),double(m2), tensor_element}));
                    }
                    }
                }
                }
            }
            }

        //prints how much was spend for running calculations
        auto end = std::chrono::high_resolution_clock::now();
        double elapsed_seconds = std::chrono::duration<double>(end - start).count();
        std::cout << "Elapsed time: " << elapsed_seconds << " seconds,\t g = " << g << "\t l = " << length << std::endl;

        //saves matrix elements in indeces-values fashion. converting output to matrices takes place in attached python file "matrix.py"
        std::ostringstream oss;
        oss << "Elapsed time: " << elapsed_seconds << "\tg = " << g << "\tl = " << length;
        std::string extrainfo = oss.str();
        functions Function;
        Function.write_to_file("output.txt", parameters, values, extrainfo);


        return 0;
    }
    
    std::cout << "n1\t" << n1 << "\tm1\t" << m1 << "\tu1\t" << u1 << std::endl;
    CompleteIntegral integral(limit,epsabs,epsrel);
    integral.changeMultiIndex(n1, m1, u1, n2, m2, u2);
    integral.integrate_norm();//sets norm for given parameters
    integral.test_integral();
    integral.print_result();

    return 0;
}

/* Example output in output.txt file located in working directory:

Elapsed time: 0.0179845	g = 1
n1	n2	u1	u2	m1	m2	integral_value
0	0	0	0	0	0	4
0	1	0	0	0	0	41.4784
0	2	0	0	0	0	159.914
1	0	0	0	0	0	41.4784
1	1	0	0	0	0	161.914
1	2	0	0	0	0	199.392
2	0	0	0	0	0	159.914
2	1	0	0	0	0	199.392
2	2	0	0	0	0	635.655

*/