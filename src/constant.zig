//! Constants used in the project from CODATA 2022 (https://physics.nist.gov/cuu/Constants/index.html).

const std = @import("std");

pub const Eh  = 4.359744722206e-18; // J
pub const a0  = 5.29177210544e-11;  // m
pub const c   = 299792458.0;        // m/s
pub const amu = 1.66053906892e-27;  // kg
pub const e   = 1.602176634e-19;    // C
pub const h   = 6.62607015e-34;     // J.s
pub const kB  = 1.380649e-23;       // J/K
pub const me  = 9.109383701528e-31; // kg
pub const NA  = 6.02214076e23;      // mol^-1

pub const A2AU = 1e-10 / a0;         // Å to a.u. of length
pub const AU2EV = Eh / e;            // a.u. of energy to eV
pub const EV2RCM = 1e-2 * e / h / c; // eV of a photon to its wavenumber in cm^-1
pub const J2AU = 1.0 / Eh;           // J to a.u. of energy
pub const U2AU = amu / me;           // g/mol to a.u

/// Atomic symbol to atomic number map.
pub const SM2AN = std.StaticStringMap(u32).initComptime(.{
    .{ "H",   1},                                                                                                                                                                                                                                 .{"He",   2},
    .{"Li",   3},                                                                                                                                             .{"Be",   4}, .{ "B",   5}, .{ "C",   6}, .{ "N",   7}, .{ "O",   8}, .{ "F",   9}, .{"Ne",  10},
    .{"Na",  11},                                                                                                                                             .{"Mg",  12}, .{"Al",  13}, .{"Si",  14}, .{ "P",  15}, .{ "S",  16}, .{"Cl",  17}, .{"Ar",  18},
    .{ "K",  19}, .{"Ca",  20}, .{"Sc",  21}, .{"Ti",  22}, .{ "V",  23}, .{"Cr",  24}, .{"Mn",  25}, .{"Fe",  26}, .{"Co",  27}, .{"Ni",  28}, .{"Cu",  29}, .{"Zn",  30}, .{"Ga",  31}, .{"Ge",  32}, .{"As",  33}, .{"Se",  34}, .{"Br",  35}, .{"Kr",  36},
    .{"Rb",  37}, .{"Sr",  38}, .{"Y",   39}, .{"Zr",  40}, .{"Nb",  41}, .{"Mo",  42}, .{"Tc",  43}, .{"Ru",  44}, .{"Rh",  45}, .{"Pd",  46}, .{"Ag",  47}, .{"Cd",  48}, .{"In",  49}, .{"Sn",  50}, .{"Sb",  51}, .{"Te",  52}, .{"I",   53}, .{"Xe",  54},
    .{"Cs",  55}, .{"Ba",  56}, .{"La",  57}, .{"Hf",  72}, .{"Ta",  73}, .{"W",   74}, .{"Re",  75}, .{"Os",  76}, .{"Ir",  77}, .{"Pt",  78}, .{"Au",  79}, .{"Hg",  80}, .{"Tl",  81}, .{"Pb",  82}, .{"Bi",  83}, .{"Po",  84}, .{"At",  85}, .{"Rn",  86},
    .{"Fr",  87}, .{"Ra",  88}, .{"Ac",  89}, .{"Rf", 104}, .{"Db", 105}, .{"Sg", 106}, .{"Bh", 107}, .{"Hs", 108}, .{"Mt", 109}, .{"Ds", 110}, .{"Rg", 111}, .{"Cn", 112}, .{"Nh", 113}, .{"Fl", 114}, .{"Mc", 115}, .{"Lv", 116}, .{"Ts", 117}, .{"Og", 118},

                                              .{"Ce",  58}, .{"Pr",  59}, .{"Nd",  60}, .{"Pm",  61}, .{"Sm",  62}, .{"Eu",  63}, .{"Gd",  64}, .{"Tb",  65}, .{"Dy",  66}, .{"Ho",  67}, .{"Er",  68}, .{"Tm",  69}, .{"Yb",  70}, .{"Lu",  71},
                                              .{"Th",  90}, .{"Pa",  91}, .{"U",   92}, .{"Np",  93}, .{"Pu",  94}, .{"Am",  95}, .{"Cm",  96}, .{"Bk",  97}, .{"Cf",  98}, .{"Es",  99}, .{"Fm", 100}, .{"Md", 101}, .{"No", 102}, .{"Lr", 103},
});

/// Atomic number to mass array.
pub const AN2M = [_]f64{
     1.007840000, // H
     4.002602000, // He
     6.938000000, // Li
     9.012183100, // Be
    10.806000000, // B
    12.009600000, // C
    14.006430000, // N
    15.999030000, // O
    18.998403163, // F
    20.179700000, // Ne
    22.989769280, // Na
    24.305500000, // Mg
    26.981538500, // Al
    28.085000000, // Si
    30.973761998, // P
    32.067000000, // S
    35.453300000, // Cl
    39.948000000, // Ar
    39.098300000, // K
    40.078000000, // Ca
    44.955908000, // Sc
    47.867000000, // Ti
    50.941500000, // V
    51.996100000, // Cr
    54.938044000, // Mn
    55.845000000, // Fe
    58.933194000, // Co
    58.693400000, // Ni
    63.546000000, // Cu
    65.380000000, // Zn
    69.723000000, // Ga
    72.630000000, // Ge
    74.921595000, // As
    78.971000000, // Se
    79.904000000, // Br
    83.798000000, // Kr
};

/// Function to revert the SM2AN map to get the symbol from the atomic number.
pub fn AN2SM(AN: u32) ![]const u8 {
    for (SM2AN.keys(), SM2AN.values()) |key, value| if (value == AN) {
        return key;
    };

    return error.InvalidAtomicNumber;
}
