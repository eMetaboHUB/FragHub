// Passerelle C++ vers RDKit pour FragHub (essai).
//
// Reproduit EXACTEMENT l'ancien rdkit_worker.process_single (scripts/rdkit_worker.py, supprimé), avec les mêmes
// fonctions RDKit appelées dans le même ordre et avec les mêmes paramètres par défaut que
// les fonctions Python correspondantes :
//   Chem.MolFromSmiles / MolFromInchi, GetNumHeavyAtoms, CalcExactMolWt, GetFormalCharge,
//   rdMolStandardize.Uncharger().uncharge, rdMolStandardize.TautomerEnumerator().Canonicalize,
//   MolToInchi, MolToInchiKey, MolToSmiles, CalcMolFormula, MolWt.
// La mise en forme des masses en texte (str(float) en Python) est faite côté Rust.

#include "fraghub_rust/src/rdkit_native/shim.h"
#include "fraghub_rust/src/rdkit_native/mod.rs.h"

#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/MolStandardize/Charge.h>
#include <GraphMol/MolStandardize/Tautomer.h>
#include <GraphMol/inchi.h>
#include <RDGeneral/RDLog.h>

#include <cmath>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace fraghub {

using RDKit::ROMol;

// --- Table des adduits : copie exacte de ADDUCT_TABLE (ancien rdkit_worker.py), mêmes expressions ---
namespace {
const double H = 1.007276, Na = 22.989218, K = 38.963158, Li = 6.015123, NH4 = 18.034164,
             H2O = 18.010565, NH3 = 17.026549, ACN = 41.026549;

struct Adduct { int mult; double delta; int charge; };

const std::map<std::string, Adduct> &adduct_table() {
  static const std::map<std::string, Adduct> t = {
      {"[M+H]+", {1, H, 1}}, {"[M+Na]+", {1, Na, 1}}, {"[M+K]+", {1, K, 1}},
      {"[M+NH4]+", {1, NH4, 1}}, {"[M+Li]+", {1, Li, 1}}, {"[M+2H]2+", {1, 2 * H, 2}},
      {"[M+H+Na]2+", {1, H + Na, 2}}, {"[M+H+K]2+", {1, H + K, 2}}, {"[M+2Na]2+", {1, 2 * Na, 2}},
      {"[M+2K]2+", {1, 2 * K, 2}}, {"[M+2Na-H]2+", {1, 2 * Na - H, 2}}, {"[M+3H]3+", {1, 3 * H, 3}},
      {"[M+2H+Na]3+", {1, 2 * H + Na, 3}}, {"[M+2H+K]3+", {1, 2 * H + K, 3}},
      {"[M+H+2Na]3+", {1, H + 2 * Na, 3}}, {"[2M+H]+", {2, H, 1}}, {"[2M+Na]+", {2, Na, 1}},
      {"[2M+K]+", {2, K, 1}}, {"[2M+NH4]+", {2, NH4, 1}}, {"[2M+2H]2+", {2, 2 * H, 2}},
      {"[2M-2H+H]+", {2, -H, 1}}, {"[2M-H2O+H]+", {2, H - H2O, 1}},
      {"[2M-2H2O+H]+", {2, H - 2 * H2O, 1}}, {"[3M+H]+", {3, H, 1}}, {"[3M+Na]+", {3, Na, 1}},
      {"[3M+K]+", {3, K, 1}}, {"[3M+NH4]+", {3, NH4, 1}}, {"[M-H2O+H]+", {1, H - H2O, 1}},
      {"[M-2H2O+H]+", {1, H - 2 * H2O, 1}}, {"[M-3H2O+H]+", {1, H - 3 * H2O, 1}},
      {"[M-4H2O+H]+", {1, H - 4 * H2O, 1}}, {"[M-5H2O+H]+", {1, H - 5 * H2O, 1}},
      {"[M+H-H2O]+", {1, H - H2O, 1}}, {"[M+H-2H2O]+", {1, H - 2 * H2O, 1}},
      {"[M+H-3H2O]+", {1, H - 3 * H2O, 1}}, {"[M-H2O+Na]+", {1, Na - H2O, 1}},
      {"[M-H2O+NH4]+", {1, NH4 - H2O, 1}}, {"[M-2H2O+NH4]+", {1, NH4 - 2 * H2O, 1}},
      {"[M-2H2O+Na]+", {1, Na - 2 * H2O, 1}}, {"[M-NH3+H]+", {1, H - NH3, 1}},
      {"[M+H-NH3]+", {1, H - NH3, 1}}, {"[M-H2+H]+", {1, H - 2 * 1.007825, 1}},
      {"[M-H2O+2H]2+", {1, 2 * H - H2O, 2}}, {"[M-3H2O+2H]2+", {1, 2 * H - 3 * H2O, 2}},
      {"[M-2H2O+2H]2+", {1, 2 * H - 2 * H2O, 2}}, {"[M+C2H3N+H]+", {1, ACN + H, 1}},
      {"[M+C2H3N+NH4]+", {1, ACN + NH4, 1}}, {"[M+C2H3N+2H]2+", {1, ACN + 2 * H, 2}},
      {"[M+3C2H3N+2H]+", {1, 3 * ACN + 2 * H, 1}}, {"[M+C2H6OS+H]+", {1, 79.021299 + H - H, 1}},
      {"[M+CH2O2+H]+", {1, 47.013305, 1}}, {"[M+CO2H2+H]+", {1, 47.013305, 1}},
      {"[M+CH4O+H]+", {1, 33.034164, 1}}, {"[M-CH4O+H]+", {1, H - 32.026215, 1}},
      {"[M]+", {1, -0.000549, 1}}, {"[M]*+", {1, -0.000549, 1}}, {"[M]2+", {1, -2 * 0.000549, 2}},
      {"[M]", {1, 0.0, 1}},
  };
  return t;
}

// calc_expected_mz : renvoie false si l'adduit est inconnu (None en Python).
// `ptype` est déjà débarrassé de ses espaces (str(precursor_type).strip(), fait côté Rust).
bool calc_expected_mz(double neutral_mass, const std::string &ptype, double &out) {
  if (ptype.empty()) return false;  // "if not precursor_type: return None"
  auto it = adduct_table().find(ptype);
  if (it == adduct_table().end()) return false;
  out = (it->second.mult * neutral_mass + it->second.delta) / it->second.charge;
  return true;
}

// Chem.MolFromSmiles : renvoie nullptr en cas d'échec (None en Python)
std::unique_ptr<ROMol> mol_from_smiles(const std::string &s) {
  try {
    RDKit::SmilesParserParams ps;
    return std::unique_ptr<ROMol>(RDKit::SmilesToMol(s, ps));
  } catch (...) {
    return nullptr;
  }
}

// Chem.MolFromInchi(inchi) = rdinchi.InchiToMol(inchi, sanitize=True, removeHs=True)
std::unique_ptr<ROMol> mol_from_inchi(const std::string &s) {
  try {
    RDKit::ExtraInchiReturnValues rv;
    return std::unique_ptr<ROMol>(RDKit::InchiToMol(s, rv, true, true));
  } catch (...) {
    return nullptr;
  }
}

// Outils créés une fois par thread (en Python : une fois par molécule, résultat identique)
struct Tools {
  RDKit::MolStandardize::Uncharger uncharger;
  RDKit::MolStandardize::TautomerEnumerator enumerator;
};
Tools &tools() {
  thread_local Tools t;
  return t;
}
}  // namespace

void init_rdkit() {
  // RDLogger.DisableLog('rdApp.*')
  boost::logging::disable_logs("rdApp.*");
}

// Descriptors.ExactMolWt(Chem.MolFromInchi(...) ou Chem.MolFromSmiles(...)), utilisé pour recalculer un m/z
// précurseur manquant (normalizer/missing_precursormz_re_calculation.rs). NaN si la molécule ne se lit pas.
double exact_mass(rust::Str mol_rs, bool is_inchi) {
  const std::string s(mol_rs.data(), mol_rs.size());
  try {
    std::unique_ptr<ROMol> mol = is_inchi ? mol_from_inchi(s) : mol_from_smiles(s);
    if (!mol) return std::nan("");
    return RDKit::Descriptors::calcExactMW(*mol);  // ExactMolWt(mol) = CalcExactMolWt(mol, onlyHeavy=False)
  } catch (...) {
    return std::nan("");
  }
}

MolResult process_single(rust::Str clean_mol_rs, rust::Str ptype_rs, double precursormz,
                         bool has_mz) {
  MolResult r{};
  const std::string clean_mol(clean_mol_rs.data(), clean_mol_rs.size());
  const std::string ptype(ptype_rs.data(), ptype_rs.size());
  try {
    std::unique_ptr<ROMol> mol = clean_mol.find("InChI=") != std::string::npos
                                     ? mol_from_inchi(clean_mol)
                                     : mol_from_smiles(clean_mol);
    if (!mol) return r;

    // --- Résolution des fragments (SMILES contenant un '.') ---
    if (clean_mol.find('.') != std::string::npos && has_mz) {
      struct Frag { std::unique_ptr<ROMol> m; unsigned n; double mass; };
      std::vector<Frag> fragments;
      size_t start = 0;
      while (true) {  // clean_mol.split('.') (garde les morceaux vides, comme Python)
        size_t dot = clean_mol.find('.', start);
        std::string part = clean_mol.substr(start, dot == std::string::npos ? std::string::npos : dot - start);
        auto m = mol_from_smiles(part);
        if (m) {
          unsigned n = m->getNumHeavyAtoms();
          double mass = RDKit::Descriptors::calcExactMW(*m);
          fragments.push_back({std::move(m), n, mass});
        }
        if (dot == std::string::npos) break;
        start = dot + 1;
      }
      if (!fragments.empty()) {
        unsigned long sum = 0, mn = fragments[0].n;
        for (auto &f : fragments) { sum += f.n; if (f.n < mn) mn = f.n; }
        if (sum > mn) {
          double mass_full = RDKit::Descriptors::calcExactMW(*mol);
          double mz_full;
          bool full_matches = calc_expected_mz(mass_full, ptype, mz_full) &&
                              std::fabs(mz_full - precursormz) <= 0.02;
          if (!full_matches) {
            std::vector<size_t> matching;
            for (size_t i = 0; i < fragments.size(); ++i) {
              double mz_f;
              if (calc_expected_mz(fragments[i].mass, ptype, mz_f) &&
                  std::fabs(mz_f - precursormz) <= 0.02)
                matching.push_back(i);
            }
            if (matching.size() == 1) {
              mol = std::move(fragments[matching[0]].m);
            } else if (matching.size() > 1) {
              // max(..., key=n) : le PREMIER fragment ayant le plus d'atomes lourds
              size_t best = matching[0];
              for (size_t i : matching) if (fragments[i].n > fragments[best].n) best = i;
              mol = std::move(fragments[best].m);
            } else {
              // sorted(..., reverse=True) est stable : à égalité, l'ordre d'origine est gardé
              std::vector<size_t> order(fragments.size());
              for (size_t i = 0; i < order.size(); ++i) order[i] = i;
              std::stable_sort(order.begin(), order.end(),
                               [&](size_t a, size_t b) { return fragments[a].n > fragments[b].n; });
              unsigned second_n = order.size() > 1 ? fragments[order[1]].n : 0;
              if (second_n <= 4) mol = std::move(fragments[order[0]].m);
            }
          }
        }
      }
    }

    // --- Neutralisation de la charge ---
    try {
      if (RDKit::MolOps::getFormalCharge(*mol) != 0) {
        std::unique_ptr<ROMol> neutral(tools().uncharger.uncharge(*mol));
        if (neutral && RDKit::MolOps::getFormalCharge(*neutral) == 0) mol = std::move(neutral);
      }
    } catch (...) {
    }

    // --- Tautomère canonique ---
    try {
      std::unique_ptr<ROMol> canon(tools().enumerator.canonicalize(*mol));
      if (canon) mol = std::move(canon);
    } catch (...) {
    }

    // --- Sorties (chacune protégée, comme les try/except du Python) ---
    try { RDKit::ExtraInchiReturnValues rv; r.inchi = rust::String(RDKit::MolToInchi(*mol, rv, "")); r.has_inchi = true; } catch (...) {}
    try { r.inchikey = rust::String(RDKit::MolToInchiKey(*mol, "")); r.has_inchikey = true; } catch (...) {}
    try { r.smiles = rust::String(RDKit::MolToSmiles(*mol)); r.has_smiles = true; } catch (...) {}
    try { r.formula = rust::String(RDKit::Descriptors::calcMolFormula(*mol)); r.has_formula = true; } catch (...) {}
    try { r.exactmass = RDKit::Descriptors::calcExactMW(*mol); r.has_exactmass = true; } catch (...) {}
    try { r.averagemass = RDKit::Descriptors::calcAMW(*mol); r.has_averagemass = true; } catch (...) {}
  } catch (...) {
  }
  return r;
}

}  // namespace fraghub
