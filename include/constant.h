#pragma once

#include   <glm/glm.hpp>
#include        <string>
#include <unordered_map>

#define OCCUPIED_ORBITALS "abcdefghijklmnopqrstuvwxyz"
#define VIRTUAL_ORBITALS  "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

#define ANGSTROM_TO_BOHR 1.889726124626

struct Atom {
    float radius, covalent; glm::vec3 color; double mass;
};

inline std::unordered_map<std::string, Atom> ptable = {
    { "H", {.radius = 053.0f, .covalent = 032.0f, .color = {255.0f / 255.0f, 255.0f / 255.0f, 255.0f / 255.0f}, .mass = 001.0078 }},
    {"He", {.radius = 031.0f, .covalent = 046.0f, .color = {217.0f / 255.0f, 255.0f / 255.0f, 255.0f / 255.0f}, .mass = 004.0026 }},
    {"Li", {.radius = 167.0f, .covalent = 133.0f, .color = {204.0f / 255.0f, 128.0f / 255.0f, 255.0f / 255.0f}, .mass = 006.9410 }},
    {"Be", {.radius = 112.0f, .covalent = 102.0f, .color = {194.0f / 255.0f, 255.0f / 255.0f, 000.0f / 255.0f}, .mass = 009.0122 }},
    { "C", {.radius = 067.0f, .covalent = 075.0f, .color = {144.0f / 255.0f, 144.0f / 255.0f, 144.0f / 255.0f}, .mass = 012.0110 }},
    { "N", {.radius = 056.0f, .covalent = 071.0f, .color = {048.0f / 255.0f, 080.0f / 255.0f, 248.0f / 255.0f}, .mass = 014.0067 }},
    { "O", {.radius = 048.0f, .covalent = 063.0f, .color = {255.0f / 255.0f, 013.0f / 255.0f, 013.0f / 255.0f}, .mass = 015.9994 }},
    { "F", {.radius = 042.0f, .covalent = 064.0f, .color = {144.0f / 255.0f, 224.0f / 255.0f, 080.0f / 255.0f}, .mass = 018.9984 }},
    {"Ne", {.radius = 038.0f, .covalent = 067.0f, .color = {179.0f / 255.0f, 227.0f / 255.0f, 245.0f / 255.0f}, .mass = 020.1797 }},
    { "P", {.radius = 098.0f, .covalent = 111.0f, .color = {255.0f / 255.0f, 128.0f / 255.0f, 000.0f / 255.0f}, .mass = 030.9738 }},
    { "S", {.radius = 088.0f, .covalent = 103.0f, .color = {255.0f / 255.0f, 255.0f / 255.0f, 048.0f / 255.0f}, .mass = 032.0650 }},
    {"Cl", {.radius = 079.0f, .covalent = 099.0f, .color = {031.0f / 255.0f, 240.0f / 255.0f, 031.0f / 255.0f}, .mass = 035.4530 }},
    {"Xe", {.radius = 108.0f, .covalent = 131.0f, .color = {066.0f / 255.0f, 158.0f / 255.0f, 176.0f / 255.0f}, .mass = 131.2930 }},
    {"El", {.radius = 020.0f, .covalent = 020.0f, .color = {066.0f / 255.0f, 158.0f / 255.0f, 176.0f / 255.0f}, .mass = 000.0000 }}
};

inline std::unordered_map<int, std::string> an2sm = {
    { 1,  "H" },
    { 2,  "He"},
    { 6,  "C" },
    { 7,  "N" },
    { 8,  "O" },
    { 9,  "F" },
    {15,  "P" },
    {16,  "S" },
    {17, "Cl" }
};

inline std::unordered_map<std::string, int> sm2an = {
    { "H",  1},
    {"He",  2},
    { "C",  6},
    { "N",  7},
    { "O",  8},
    { "F",  9},
    { "P", 15},
    { "S", 16},
    {"Cl", 17}
};
