#pragma once

#define A2BOHR 1.889726124626
#define BOHR2A 0.529177210903
#define AMU2AU 1822.888486192

#include "argparse.h"

inline std::unordered_map<int, std::string> AN2SM = {
    { 1,  "H"},
    { 2, "He"},
    { 6,  "C"},
    { 7,  "N"},
    { 8,  "O"},
    { 9,  "F"},
    {10, "Ne"},
    {11, "Na"},
    {14, "Si"},
    {15,  "P"},
    {16,  "S"},
    {17, "Cl"},
    {18, "Ar"}
};

inline std::unordered_map<std::string, int> SM2AN = {
    { "H",  1,},
    {"He",  2,},
    { "C",  6,},
    { "N",  7,},
    { "O",  8,},
    { "F",  9,},
    {"Ne", 10,},
    {"Na", 11,},
    {"Si", 14,},
    { "P", 15,},
    { "S", 16,},
    {"Cl", 17,},
    {"Ar", 18,}
};
