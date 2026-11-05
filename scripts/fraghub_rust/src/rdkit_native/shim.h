// En-tête de la passerelle C++ vers RDKit (voir shim.cc)
#pragma once
#include "rust/cxx.h"

namespace fraghub {
struct MolResult;
void init_rdkit();
MolResult process_single(rust::Str clean_mol, rust::Str ptype, double precursormz, bool has_mz);
double exact_mass(rust::Str mol, bool is_inchi);
}  // namespace fraghub
