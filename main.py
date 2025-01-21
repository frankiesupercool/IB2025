import os
import subprocess
from rpy2.robjects import r, pandas2ri
import rpy2.robjects.packages as rpackages
import rpy2.robjects as robjects

# Enable the conversion between pandas DataFrames and R data frames
pandas2ri.activate()

# TODO Rmpi maybe not needed anymore
# install.packages(c("varSelRF", "data.table", "ggplot2", "dplyr", "foreach", "snow", "parallel", "Rmpi", "randomForest))
#
# # Install Bioconductor and required packages
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#     install.packages("BiocManager")
# }
# BiocManager::install(c("RTCGA", "RTCGA.mRNA", "RTCGA.methylation", "limma"))


# TODO python script not usable atm - use setup above in R
def install_r_packages(packages):
    """Install required R packages."""
    utils = rpackages.importr("utils")
    utils.chooseCRANmirror(ind=1)  # Select the first CRAN mirror
    for package in packages:
        if not rpackages.isinstalled(package):
            utils.install_packages(package)


def setup_environment():
    """Ensure R dependencies are installed."""
    r_dependencies = [
        "varSelRF",
        "data.table",
        "ggplot2",
        "dplyr",
        "foreach",
        "snow",
        "parallel",
        "Rmpi",
        "randomForest",
        "BiocManager"
    ]
    install_r_packages(r_dependencies)
    # Install Bioconductor packages if needed
    r('''
    if (!requireNamespace("RTCGA", quietly = TRUE)) {
        BiocManager::install("RTCGA", "RTCGA.mRNA", "RTCGA.methylation", "limma")
    }
    ''')


def run_r_script(script_name):
    """Run an R script."""
    try:
        print(f"Running R script: {script_name}")
        subprocess.check_call(['Rscript', script_name])
    except subprocess.CalledProcessError as e:
        print(f"Error while running {script_name}: {e}")


def main():
    # Step 1: Setup environment
    print("Setting up R environment...")
    setup_environment()

    # Step 2: Run the R scripts in sequence
    r_scripts = [
        "download_tcga_data.R",
        "process_tcga_data.R",
        "remove.false.normal.labels.R",
        "merge.normal.samples.R",
        "varSelRF.R"
    ]

    for script in r_scripts:
        if os.path.exists(script):
            run_r_script(script)
        else:
            print(f"Script {script} not found. Ensure all scripts are in the current directory.")

    print("All scripts have been executed.")


if __name__ == "__main__":
    main()