#ifndef Flag_parser
#define Flag_parser

#include <vector>

class FlagParser{
    public:
        int argc;
        char** argv;

        bool debugFlag;
        std::vector<int> debugValues;
        
        bool change_gFlag;
        double change_gValue;

        bool minmaxFlag;
        double nmin, nmax, mmin, mmax, umin, umax;

        bool change_lFlag;
        double change_lValue;

        FlagParser(int argc, char** argv){
            this->argc = argc;
            this->argv = argv;
            this->debugFlag = false;
            this->change_gFlag = false;
            this->minmaxFlag = false;
            this->change_lFlag = false;
        }
        ~FlagParser(){}
        void parse_flags();
    };

#endif