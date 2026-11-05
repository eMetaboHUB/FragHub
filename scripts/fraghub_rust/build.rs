// Compile la passerelle C++ vers RDKit (src/rdkit_native/shim.cc) et la lie aux bibliothèques
// RDKit précompilées de conda-forge (paquet librdkit-dev + libboost-headers).
//
// Dossier de l'environnement conda : variable RDKIT_PREFIX, sinon ~/rdkit_cpp
// (ex. sous Windows : mamba create -p C:\Users\<nom>\rdkit_cpp -c conda-forge librdkit-dev libboost-headers).
use std::path::PathBuf;

fn main() {
    println!("cargo:rerun-if-changed=src/rdkit_native/shim.cc");
    println!("cargo:rerun-if-changed=src/rdkit_native/shim.h");
    println!("cargo:rerun-if-env-changed=RDKIT_PREFIX");

    let prefix = match std::env::var("RDKIT_PREFIX") {
        Ok(p) if !p.is_empty() => PathBuf::from(p),
        _ => {
            let home = std::env::var("USERPROFILE").or_else(|_| std::env::var("HOME")).unwrap_or_default();
            PathBuf::from(home).join("rdkit_cpp")
        }
    };
    let windows = std::env::var("CARGO_CFG_TARGET_OS").unwrap() == "windows";
    // Sous Windows, conda range en-têtes et bibliothèques dans Library\
    let base = if windows { prefix.join("Library") } else { prefix.clone() };
    let include = base.join("include");
    if !include.join("rdkit").join("GraphMol").exists() {
        panic!(
            "RDKit C++ introuvable dans {}.\n\
             Installez-le avec : mamba create -p {} -c conda-forge librdkit-dev=2026.09.1 libboost-headers\n\
             ou indiquez son dossier dans la variable RDKIT_PREFIX.",
            prefix.display(), prefix.display()
        );
    }

    let mut b = cxx_build::bridge("src/rdkit_native/mod.rs");
    b.file("src/rdkit_native/shim.cc")
        .include(&include)
        .include(include.join("rdkit"))
        .std("c++20")
        .warnings(false);
    if windows {
        // RDKIT_DYN_LINK : importer les fonctions depuis les DLL de RDKit (__declspec(dllimport))
        // BOOST_ALL_NO_LIB : pas de liaison automatique à Boost (les DLL de RDKit chargent les leurs)
        // /Zc:__cplusplus : MSVC annonce sa vraie version du C++
        b.flag("/EHsc").flag("/bigobj").flag("/Zc:__cplusplus").flag("/permissive-")
            .define("NOMINMAX", None).define("RDKIT_DYN_LINK", None).define("BOOST_ALL_NO_LIB", None);
    }
    b.compile("fraghub_rdkit_shim");

    let libdir = base.join("lib");
    println!("cargo:rustc-link-search=native={}", libdir.display());
    for lib in ["RDKitMolStandardize", "RDKitRDInchiLib", "RDKitDescriptors", "RDKitSmilesParse",
                "RDKitSubstructMatch", "RDKitGraphMol", "RDKitRDGeneral"] {
        println!("cargo:rustc-link-lib=dylib={lib}");
    }
    if !windows {
        // Linux / macOS : chemin des bibliothèques gravé dans le module (pour le développement)
        println!("cargo:rustc-link-arg=-Wl,-rpath,{}", libdir.display());
    }
}
