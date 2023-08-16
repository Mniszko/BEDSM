#include <iostream>
#include <string>
#include <vector>
#include <sstream>

#include "flag_parser.h"


void FlagParser::parse_flags(){
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--debug") {
            this->debugFlag = true;
            if (i + 3 <= argc) {
                for (int j = 0; j < 3; ++j) {
                    this->debugValues.push_back(std::stoi(argv[i + j + 1]));
                }
                i += 3;
            } else if (i + 6 <= argc) {
                for (int j = 0; j < 6; ++j) {
                    this->debugValues.push_back(std::stoi(argv[i + j + 1]));
                }
                i += 6;
            } else {
                std::cerr << "Error: --debug expects 3 or 6 integer arguments." << std::endl;
                return;
            }
        } else if (arg == "--change_g" && i + 1 < argc) {
            this->change_gFlag = true;
            this->change_gValue = std::stod(argv[i + 1]);
            i += 1;
        }else if (arg == "--change_l" && i + 1 < argc) {
            this->change_lFlag = true;
            this->change_lValue = std::stod(argv[i + 1]);
            i += 1;
        } else if (arg == "--minmax" && i + 6 < argc) {
            this->minmaxFlag = true;
            this->nmin = std::stod(argv[i + 1]);
            this->nmax = std::stod(argv[i + 2]);
            this->mmin = std::stod(argv[i + 3]);
            this->mmax = std::stod(argv[i + 4]);
            this->umin = std::stod(argv[i + 5]);
            this->umax = std::stod(argv[i + 6]);
            i += 6;
        } else {
            std::cerr << "Unknown or incorrectly used flag." << std::endl;
            return;
        }
    }
}