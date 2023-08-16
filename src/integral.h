#ifndef INTEGRAL_H
#define INTEGRAL_H

class CompleteIntegral{
    private:
        int m1, n1, u1, n2, m2, u2, limit;
        double result, error, epsabs, epsrel;
        double length, rev_length;
        double norm_psi(double, void*);
        double norm1, norm2;
        int key;
        static void check_if_valid_key(int k);
    public:
        CompleteIntegral(int limit, double epsabs, double epsrel){
            this->epsabs = epsabs;
            this->epsrel = epsrel;
            this->limit = limit;
            this->m1 = 0;
            this->n1 = 0;
            this->u1 = 0;
            this->m2 = 0;
            this->n2 = 0;
            this->u2 = 0;
            this->result = 0;
            this->error = 0;
            this->norm1 = 0;
            this->norm2 = 0;
            this->length = 1;
            this->rev_length = 1/this->length;
            this->key = 1;
            this->check_if_valid_key(this->key);
        }
        void change_key(int k){
            this->key = k;
            this->check_if_valid_key(this->key);
        }
        void add_to_key(){
            this->key += 1;
            this->check_if_valid_key(this->key);
        }
        void change_length(double l){
            this->length = l;
            this-> rev_length = 1/l;
        }
        ~CompleteIntegral(){}
        double get_result(){return this->result;}
        void integrate_norm();
        double integrate_norm(double, double, double);
        void changeMultiIndex(int, int, int, int, int, int);//m should be always 0 or positive (negative one isn't implemented yet, it doesn't require much time, but can lead to unexpected errors), n and u can be either positive or negative or zero 
        void changeLimit(int);
        double give_norm(int n){
            if (n==n1){
                return norm1;
            } else if (n==n2){
                return norm2;
            }
            return 0;
        } //helper functions that shouldn't exist, but rewriting the other ones would be time consuming
        double give_norm(int m, int u){
            if (m==m1 && u==u1){
                return norm1;
            } else if (m==m2 && u==u2){
                return norm2;
            }
            return 0;
        }
        void restart_counter(){
            this->result = 0;
            this->error = 0;
        }
        void print_result();

        double rhoPhiLaplacianIntegral_real(int, int);
        double rhoPhiLaplacianIntegral_imaginary(int, int);
        double rhoPhiHarmonicPotentialIntegral(int, int);
        double zLaplacianIntegral_real(int, double);
        double zLaplacianIntegral_imaginary(int, double);
        double rhoPhiLaplacianIntegralWithHarmonic_real(int m, int u);
        double rhoPhiLaplacianIntegralWithHarmonic_imaginary(int m, int u);
        double fastLaplacianWithHarmonicEnergy();

        double fast_add_over_harmonic();
        double test_integral();
        double Integrate_over_harmonic();
        double Integrate_over_harmonic_with_delta();//3D movement case without magnetic moment interactions
        double Integrate_over_delta_1D_movement_real(double);
        double Integrate_over_delta_1D_movement_imaginary(double);
        double Integrate_over_harmonic_with_delta_potential_1D_motion_head_crush(double);
        double* Integrate_over_harmonic_with_delta_1D_movement(double);//replication of Lieb-Lienieger energy result in relation to variable g. All quantum numbers should be the same in this case, except n, for which n1 = -n2 or n1 != n2. First case is quite easy

};


#endif /* INTEGRAL_H */