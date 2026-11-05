# build_dataset

Fabrique un **petit jeu de données de test** à partir d'un dossier d'entrée FragHub, pour lancer
le pipeline en quelques minutes au lieu de plusieurs heures.

Pour chaque fichier `.msp`, `.mgf`, `.json` et `.csv` du dossier, l'outil garde un échantillon
aléatoire de spectres (tirage uniforme sur tout le fichier) et l'écrit dans un fichier de **même
nom et même format** dans le dossier de sortie.

- Les spectres sont recopiés **octet pour octet**, sans réécriture.
- Le découpage suit les mêmes règles que les lecteurs de FragHub
  (`scripts/fraghub_rust/src/convertors/loaders.rs`) :
  - MSP : blocs de lignes séparés par une ligne vide ;
  - MGF : jusqu'à `END IONS` ;
  - JSON : objets `{...}` de premier niveau (tableau GNPS ou JSONL) ;
  - CSV : en-tête + un enregistrement par spectre.
- Les spectres gardés restent dans leur ordre d'origine.
- Le tirage est reproductible : même `--seed`, même échantillon.
- Le dossier d'entrée n'est jamais modifié. L'outil refuse d'écrire à l'intérieur.

## Utilisation

Depuis la racine du projet :

```powershell
cargo run --release --manifest-path scripts/build_dataset/Cargo.toml -- "C:\chemin\input" [dossier_sortie] [--percent 20 | --per-file 2000] [--seed 42]
```

Par défaut, la sortie est `<dossier_entrée>-sample`, à côté du dossier d'entrée, avec 2000 spectres
par fichier. Un fichier `_sample_info.txt` (ignoré par FragHub) résume ce qui a été gardé.

L'outil n'a aucune dépendance : il compile en quelques secondes.
