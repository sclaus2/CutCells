#!/usr/bin/env bash
# CutCells documentation build helper.

set -e

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'
DOCS_CONDA_ENV="${CUTCELLS_DOCS_CONDA_ENV:-}"
export LC_ALL=C
export LANG=C

print_status() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

if [ ! -f "conf.py" ]; then
    print_error "This script must be run from the CutCells/docs directory"
    exit 1
fi

initialize_conda() {
    if [ -z "${CONDA_EXE:-}" ]; then
        if command -v conda > /dev/null 2>&1; then
            CONDA_EXE=$(which conda)
        elif [ -f "$HOME/miniconda3/bin/conda" ]; then
            CONDA_EXE="$HOME/miniconda3/bin/conda"
        elif [ -f "$HOME/anaconda3/bin/conda" ]; then
            CONDA_EXE="$HOME/anaconda3/bin/conda"
        else
            print_error "Conda not found. Activate an environment manually or install conda."
            exit 1
        fi
    fi

    eval "$($CONDA_EXE shell.bash hook)"
}

activate_env() {
    if [ -z "$DOCS_CONDA_ENV" ]; then
        print_status "Using the current shell environment"
        return
    fi

    print_status "Setting up conda environment: ${DOCS_CONDA_ENV}"
    initialize_conda
    conda activate "${DOCS_CONDA_ENV}"
}

install_deps() {
    activate_env
    print_status "Installing documentation dependencies..."
    python -m pip install -r requirements.txt
    print_success "Dependencies installed"
}

clean_build() {
    print_status "Cleaning documentation build directory..."
    rm -rf _build/
    print_success "Build directory cleaned"
}

build_docs() {
    activate_env
    print_status "Building HTML documentation..."
    sphinx-build -b html . _build/html
    print_success "Documentation built at _build/html/index.html"
}

serve_docs() {
    local port=${1:-8000}
    if [ ! -f "_build/html/index.html" ]; then
        print_warning "Documentation is not built yet; building first."
        build_docs
    fi

    print_status "Serving documentation at http://localhost:${port}"
    cd _build/html
    python -m http.server "$port"
}

live_docs() {
    activate_env
    if ! command -v sphinx-autobuild > /dev/null 2>&1; then
        print_status "Installing sphinx-autobuild..."
        python -m pip install sphinx-autobuild
    fi
    sphinx-autobuild . _build/html --host 0.0.0.0 --port 8000
}

check_links() {
    activate_env
    print_status "Checking links..."
    sphinx-build -b linkcheck . _build/linkcheck
    print_success "Link check completed"
}

show_help() {
    echo "CutCells Documentation Build Script"
    echo
    echo "Usage: $0 [COMMAND]"
    echo
    echo "Commands:"
    echo "  install     Install documentation dependencies"
    echo "  clean       Clean the documentation build directory"
    echo "  build       Build HTML documentation"
    echo "  serve       Serve built documentation locally"
    echo "  live        Start a live-reload documentation server"
    echo "  check       Check links"
    echo "  all         Clean, build, and serve"
    echo "  help        Show this help message"
    echo
    echo "Environment:"
    echo "  Set CUTCELLS_DOCS_CONDA_ENV to choose a conda environment."
}

case "${1:-help}" in
    install)
        install_deps
        ;;
    clean)
        clean_build
        ;;
    build)
        build_docs
        ;;
    serve)
        port=8000
        if [ "${2:-}" = "--port" ] && [ -n "${3:-}" ]; then
            port=$3
        fi
        activate_env
        serve_docs "$port"
        ;;
    live)
        live_docs
        ;;
    check)
        check_links
        ;;
    all)
        port=8000
        if [ "${2:-}" = "--port" ] && [ -n "${3:-}" ]; then
            port=$3
        fi
        clean_build
        build_docs
        serve_docs "$port"
        ;;
    help|--help|-h)
        show_help
        ;;
    *)
        print_error "Unknown command: $1"
        echo
        show_help
        exit 1
        ;;
esac
