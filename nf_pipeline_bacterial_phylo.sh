#!/bin/bash

# Simplified script to run bacterial phylogenetic pipeline using SNP data

# required to run .nf script + "modules" should be a subdirectory
projectDir="/home/michall/git/pzh-phylogenetic-pipeline" # location of main project

# docker images required to execute this pipeline
main_image="pzh_pipeline_phylogenetic:latest" # main image used by phylogenetic pipeline
prokka_image="staphb/prokka:latest" # name of the image with prokka software

## Nextflow executor (local or slurm)
profile="local"

# location of input data - USER-provided no defaults
metadata=""
inputDir=""
inputType=""
genus="" # analyzed genus only Salmonella Escherichi and Campylobacter are currently supported

# output - localization of output + prefix added to all results
# as we aggregate multiple files we cannot "guess" it as e.g. we do for NGS pipeline
results_dir="./results"
results_prefix=""

# Pipeline-specific parameteres with defaults
# Defaults are hardcoded in this scrupt NOT in the .nf file 
model="GTR+G" # Model for raxml
startingTrees=10 # Number of random initial trees
bootstrap=200 # Number of bootstraps
minSupport=70 # Minimum support for a branch to keep it in a tree
clockRate="" # User can still overrride any built-in and estimated values fron the alignment. If empty data derived
threads=36 # ilosc watkow uzywanych maksymalnie przez pipeline
# QC params
thresholdN=100 # maksymalna ilosc N w genomie
thresholdAmbigous=100 # maksymalna ilosc znakow ambigous w genomie


# Usage function to display help
usage() {
    echo "Użycie: $0 --input_dir ŚCIEŻKA --input_type TYP --results_dir ŚCIEŻKA --profile TYP --metadata metadata.txt --genus Salmonella --results_prefix prefiks"
    echo
    echo "Skrypt uruchamia bakteryjny pipeline filogenetyczny oparty na danych SNP."
    echo
    echo "Wymagane parametry:"
    echo "  -m, --metadata ŚCIEŻKA            Ścieżka do pliku metadata (WYMAGANE)"
    echo "  -i, --inputDir ŚCIEŻKA            Ścieżka do katalogu z danymi wejściowymi SNP (WYMAGANE)"
    echo "  -t, --inputType TYP               Typ danych wejściowych: gff lub fasta (WYMAGANE)"
    echo "  -g, --genus NAZWA                 Rodzaj bakterii: Salmonella, Escherichia lub Campylobacter (WYMAGANE)"
    echo "  -p, --results_prefix PREFIKS      Prefiks dodawany do wszystkich plików wynikowych (WYMAGANE)"
    echo "  -d, --projectDir ŚCIEŻKA          Ścieżka do lokalnej wersji repozytorium"
    echo
    echo "Opcjonalne parametry:"
    echo "  -o, --results_dir ŚCIEŻKA         Katalog wynikowy (domyślnie: ./results)"
    echo "  -x, --profile NAZWA               Profil wykonania Nextflow (dozwolone: 'local' lub 'slurm', domyślnie: local)"
    echo "  -r, --clockRate WARTOŚĆ           Stała szybkość zegara molekularnego do TimeTree; jeśli nie podana, szacowana z danych lub ustawiana na wartość domyślną dla rodzaju"
    echo "  --model MODEL                     Model substytucji dla RAxML (domyślnie: GTR+G)"
    echo "  --startingTrees LICZBA            Liczba drzew startowych dla RAxML (domyślnie: 10)"
    echo "  --bootstrap LICZBA                Liczba replikacji bootstrap dla RAxML (domyślnie: 200)"
    echo "  --minSupport LICZBA               Minimalne wsparcie gałęzi, by pozostała w końcowym drzewie (domyślnie: 70)"
    echo "  --threads LICZBA                  Liczba rdzeni CPU do wykorzystania (domyślnie: 40)"
    echo "  --thresholdN LICZBA               Maksymalna liczba znaków 'N' w genomie (domyślnie: 100)"
    echo "  --thresholdAmbigous LICZBA        Maksymalna liczba niejednoznacznych znaków w genomie (domyślnie: 100)"
    echo "  --main_image NAZWA:TAG            Obraz Docker zawierający narzędzia używane przez pipeline"
    echo "  --prokka_image NAZWA:TAG          Obraz Docker z oprogramowaniem Prokka"
    echo "  -h, --help                      Show this help message"
    exit 1
}

# Parse arguments using GNU getopt
TEMP=$(getopt -o hm:i:t:g:p:d:o:x:r: \
--long metadata:,inputDir:,inputType:,genus:,results_prefix:,projectDir:,results_dir:,profile:,clockRate:,model:,startingTrees:,bootstrap:,minSupport:,threads:,thresholdN:,thresholdAmbigous:,main_image:,prokka_image:,help \
-n "$0" -- "$@")

if [ $? != 0 ]; then usage; fi
eval set -- "$TEMP"

# Parse args into variables
while true; do
  case "$1" in
    -m|--metadata) metadata="$2"; shift 2 ;;
    -i|--inputDir) inputDir="$2"; shift 2 ;;
    -t|--inputType) inputType="$2"; shift 2 ;;
    -g|--genus) genus="$2"; shift 2 ;;
    -p|--results_prefix) results_prefix="$2"; shift 2 ;;
    -d|--projectDir) projectDir="$2"; shift 2 ;;
    -o|--results_dir) results_dir="$2"; shift 2 ;;
    -x|--profile) profile="$2"; shift 2 ;;
    -r|--clockRate) clockRate="$2"; shift 2 ;;
    --model) model="$2"; shift 2 ;;
    --startingTrees) startingTrees="$2"; shift 2 ;;
    --bootstrap) bootstrap="$2"; shift 2 ;;
    --minSupport) minSupport="$2"; shift 2 ;;
    --threads) threads="$2"; shift 2 ;;
    --thresholdN) thresholdN="$2"; shift 2 ;;
    --thresholdAmbigous) thresholdAmbigous="$2"; shift 2 ;;
    --main_image) main_image="$2"; shift 2 ;;
    --prokka_image) prokka_image="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    --) shift; break ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

# Check required arguments
if [[ -z "$metadata" || -z "$inputDir" || -z "$inputType" || -z "$genus" || -z "$results_prefix" ]]; then
    echo "Błąd: brakuje wymaganych argumentów!"
    usage
fi

# 1. Check if input paths exist
if [ ! -f "$metadata" ]; then
    echo "Błąd: plik metadata '$metadata' nie istnieje."; exit 1
fi
if [ ! -d "$inputDir" ]; then
    echo "Błąd: katalog wejściowy '$inputDir' nie istnieje."; exit 1
fi
if [ ! -d "$projectDir" ]; then
    echo "Błąd: katalog projektu '$projectDir' nie istnieje."; exit 1
fi

# 2. Check if Docker images exist locally
if ! docker image inspect "$main_image" > /dev/null 2>&1; then
    echo "Błąd: obraz Docker '$main_image' nie istnieje lokalnie."; exit 1
fi
if ! docker image inspect "$prokka_image" > /dev/null 2>&1; then
    echo "Błąd: obraz Docker '$prokka_image' nie istnieje lokalnie."; exit 1
fi

# 3. Validate profile
if [[ "$profile" != "local" && "$profile" != "slurm" ]]; then
    echo "Błąd: nieprawidłowy profil Nextflow: '$profile'. Dozwolone: 'local', 'slurm'."; exit 1
fi

# 4. Validate genus
if [[ "$genus" != "Salmonella" && "$genus" != "Escherichia" && "$genus" != "Campylobacter" ]]; then
    echo "Błąd: nieprawidłowy rodzaj bakterii: '$genus'. Dozwolone: Salmonella, Escherichia, Campylobacter."; exit 1
fi

# 5. Validate inputType
if [[ "$inputType" != "fasta" && "$inputType" != "gff" ]]; then
    echo "Błąd: nieprawidłowy typ danych wejściowych: '$inputType'. Dozwolone: fasta, gff."; exit 1
fi

# 6. Validate thread count
max_cpus=$(nproc)
if ! [[ "$threads" =~ ^[0-9]+$ ]] || [ "$threads" -le 0 ] || [ "$threads" -gt "$max_cpus" ]; then
    echo "Błąd: liczba wątków '$threads' jest nieprawidłowa. Dozwolone wartości: od 1 do $max_cpus."; exit 1
fi

if [ $((threads % 12)) -ne 0 ]; then
    echo "Błąd: liczba wątków ($threads) musi być wielokrotnością 12"; exit 1
fi

nextflow run ${projectDir}/nf_pipeline_bacterial_phylo.nf \
	     --input_dir ${inputDir} \
             --input_type ${inputType} \
	     --metadata ${metadata}  \
	     --genus ${genus} \
             --input_prefix ${results_prefix} \
             --clockRate "${clockRate}" \
	     --model ${model} \
	     --starting_trees ${startingTrees} \
	     --bootstrap ${bootstrap} \
	     --min_support ${minSupport} \
	     --threshold_Ns ${thresholdN} \
	     --threshold_ambiguities ${thresholdAmbigous} \
	     --main_image ${main_image} \
	     --prokka_image ${prokka_image} \
	     --results_dir ${results_dir} \
	     --threads ${threads} \
	     -profile ${profile} \
	     -with-trace

