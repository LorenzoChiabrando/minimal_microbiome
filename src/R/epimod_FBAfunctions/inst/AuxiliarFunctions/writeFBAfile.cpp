#include <Rcpp.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>

using namespace Rcpp;
using namespace std;

/**
 * @brief Writes optimization model data to a text file.
 * 
 * This function exports a stoichiometry matrix and associated vectors into a text file,
 * formatted specifically for use by optimization software. The file includes reaction identifiers,
 * the stoichiometry matrix, objective coefficients, variable bounds, and additional metadata
 * such as maximum and average biomass estimates for the species being modeled.
 *
 * @param S NumericMatrix, the stoichiometry matrix of the system, where rows represent constraints
 *          and columns represent variables.
 * @param react_id CharacterVector, the identifiers for the reactions corresponding to the columns of S.
 * @param obj_coef NumericVector, the coefficients for the objective function.
 * @param lowbnd NumericVector, the lower bounds for the variables.
 * @param uppbnd NumericVector, the upper bounds for the variables.
 * @param rb NumericVector, the right-hand side coefficients of the constraints, assuming all are equalities.
 * @param gene_assoc NumericVector, a 0/1 flag indicating if each reaction is gene-associated (same order as react_id).
 * @param model_name string, the base name for the output file, which will be saved as <model_name>.txt.
 * @param write bool, controls whether the file is actually written; if false, the function returns immediately.
 * @param wd string, the working directory where the file will be saved.
 * @param bioMax double, the maximum biomass this model can achieve (default -1 if unspecified).
 * @param bioMean double, the average biomass for this species (default -1 if unspecified).
 * @param bioMin double, the minimum biomass for this species (default -1 if unspecified).
 *
 * @details
 * The output format is as follows (with one extra line for gene_assoc):
 * - Line 1: Space-separated reaction IDs.
 * - Line 2: Space-separated 0/1 flags for whether each reaction is gene-associated.
 * - Line 3: "nrow ; ncol ; GLP_MAX"
 * - Line 4: The total number of non-zero entries in the matrix S (nnz).
 * - Line 5: Biomax value.
 * - Line 6: Biomean value.
 * - Line 7: Biomin value.
 * - Line 8: The objective coefficients (obj_coef).
 * - Next m lines (row bounds), Next n lines (column bounds),
 * - Finally the non-zero matrix entries of S, formatted "row_index ; col_index ; value".
 */

// [[Rcpp::export]]
void writeModelCpp(NumericMatrix S,
                   CharacterVector react_id,
                   NumericVector obj_coef,
                   NumericVector lowbnd,
                   NumericVector uppbnd,
                   NumericVector rb,
                   NumericVector gene_assoc, // <--- Nuovo argomento
                   string model_name,
                   bool write,
                   string wd,
                   double bioMax = -1,
                   double bioMean = -1,
                   double bioMin = -1,
                   int pFBAFlag  = -1)
{
    if (!write) {
      Rcout << "Skipping file writing (write=false)" << std::endl;
      return;
    }

    int nrow = S.nrow();
    int ncol = S.ncol();

    // Prepara il path e apri il file
    string fileName = wd + model_name + ".txt";
    ofstream file(fileName);
    if (!file.is_open()) {
        Rcout << "Failed to open file: " << fileName << endl;
        return;
    }

    // Buffer per tutto il contenuto
    stringstream buffer;

    // 1) Scriviamo i nomi delle reazioni (unica riga, spazio-separati)
    for (int i = 0; i < react_id.size(); i++) {
        buffer << as<string>(react_id[i]);
        if (i < react_id.size() - 1) buffer << " ";
    }
    buffer << "\n";

    // 2) Scriviamo la riga con gene_assoc (unica riga, spazio-separati)
    for (int i = 0; i < gene_assoc.size(); i++) {
        buffer << gene_assoc[i];
        if (i < gene_assoc.size() - 1) buffer << " ";
    }
    buffer << "\n";

    // 3) Dimensioni e tipo di problema
    buffer << nrow << " ; " << ncol << " ; GLP_MAX" << " ; " << pFBAFlag << "\n";
    

    // 4) Prepariamo i dati della matrice S (solo i non-zero)
    int nnz = 0;
    stringstream matrixData;
    matrixData << scientific; // notazione scientifica

    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
            if (S(i, j) != 0.0) {
                matrixData << (i + 1) << " ; "
                           << (j + 1) << " ; "
                           << S(i, j) << "\n";
                nnz++;
            }
        }
    }

    // 5) Scriviamo il numero di elementi non-zero
    buffer << nnz << "\n";

    // 6) BioMax, BioMean, BioMin
    buffer << bioMax << "\n";
    buffer << bioMean << "\n";
    buffer << bioMin << "\n";

    // 7) Coefficienti obiettivo (unica riga)
    for (int i = 0; i < obj_coef.size(); i++) {
        buffer << obj_coef[i];
        if (i < obj_coef.size() - 1) buffer << " ";
    }
    buffer << "\n";

    // 8) Row bounds (constraints, dimensione nrow)
    //    ipotizziamo "GLP_FX" e i = row
    //    rb[i] come valore di eq
    for (int i = 0; i < nrow; i++) {
        buffer << "GLP_FX ; " << rb[i] << " ; " << rb[i] << "\n";
    }

    // 9) Column bounds (variables, dimensione ncol)
    //    lowbnd[j] e uppbnd[j]
    for (int j = 0; j < ncol; j++) {
        buffer << "GLP_DB ; " << lowbnd[j] << " ; " << uppbnd[j] << "\n";
    }

    // 10) Aggiungiamo le righe della matrice S non-zero
    buffer << matrixData.str();

    // Fine: scriviamo il buffer nel file
    file << buffer.str();
    file.close();
    Rcout << "File written to: " << fileName << endl;
}

