#' @title A Granularity-Based Approach to identify Spatially Variable Genes
#' @description  This function is designed to identify spatially variable genes through a granularity-based approach. 
#' @details This function utilizes a MxD matrix (Coords) representing D-dimensional coordinates with M spots and a sparse, NxM expression matrix (ExpMat_Sp) with N genes and M spots. 
#' @usage GrabSVG(Coords, ExpMat_Sp, D_1 = 1.0, D_2 = 3.0, 
#' Exp_Norm = TRUE, Coords_Norm_Method = c("Sliced", "Overall", "None"))
#' @param Coords A M x D matrix representing D-dimensional coordinates for M spots
#' @param ExpMat_Sp A sparse, N x M expression matrix in dgCMatrix class with N genes and M spots
#' @param D_1 Size of the small patch
#' @param D_2 Size of the big patch
#' @param Exp_Norm A Boolean value indicating whether the expression matrix should be normalized
#' @param Coords_Norm_Method Normalization method for the coordinates matrix, which can be "None", "Sliced", or "Overall".
#' @return A data frame with the name of genes and corresponding p-values.
#' @examples
#' Coords <- expand.grid(1:100,1:100, 1:3)
#' RandFunc <- function(n) floor(10 * stats::rbeta(n, 1, 5))
#' Raw_Exp <- Matrix::rsparsematrix(nrow = 10^4, ncol = 3*10^4, density = 0.0001, rand.x = RandFunc)
#' Filtered_ExpMat <- SpFilter(Raw_Exp)
#' rownames(Filtered_ExpMat) <- paste0("Gene_", 1:nrow(Filtered_ExpMat))
#' P_values <- GrabSVG(Coords, Filtered_ExpMat)
#' @importFrom stats quantile plnorm
#' @export

GrabSVG <- function(Coords, ExpMat_Sp, D_1 = 1.0, D_2 = 3.0, Exp_Norm = TRUE, Coords_Norm_Method = c("Sliced", "Overall", "None")){
  # define a function to calculate variance of local means
  VarLocalMeans <- function(DK){
    KDBinary <- RANN::nn2(Coords, k = 100, searchtype = "radius", radius = DK)$nn.idx
    KDBinary <- Matrix::Matrix(KDBinary, sparse = TRUE)
    KDBinary_rowp <- sparseMatrixStats::rowCounts(KDBinary>0)
    KDBinary_rowind <- rep(1:length(KDBinary_rowp), KDBinary_rowp)
    KDBinary <- spam::as.spam.dgCMatrix(KDBinary)
    PatchesCells <- Matrix::sparseMatrix(x = rep(1, length(KDBinary@entries)),
                                        i = KDBinary_rowind,
                                        j = KDBinary@entries)
    PatchesCells_Centroid <- Matrix::Diagonal(x = (Matrix::colSums(PatchesCells) > 1))
    PatchesCells <- PatchesCells - PatchesCells_Centroid
    diagMatrix_sparse <- Matrix::Diagonal(x = 1/Matrix::colSums(PatchesCells))
    X_kj_Inter <- ExpMat_SpNorm %*% PatchesCells
    X_kj <- X_kj_Inter %*% diagMatrix_sparse
    var_X_K <- sparseMatrixStats::rowVars(X_kj)
    return(var_X_K)
  }
  
  # normalization of expressions
  SpMinMax <- function(ExpMat_Sp){
    # will use approximation
    RowMaxs <- sparseMatrixStats::rowMaxs(ExpMat_Sp)
    Norm_MaxMat <- Matrix::Diagonal(x = 1/RowMaxs)
    ExpMat_Norm <- Norm_MaxMat %*% ExpMat_Sp
  }
  
  # Estimate log-norm parameters
  Lnorm_estimation <- function(TSMatrix){
    T_matrix_sum_upper90 <- quantile(TSMatrix, 0.90)
    T_matrix_sum_mid <- TSMatrix[which(TSMatrix<T_matrix_sum_upper90)]
    LnormPar <- fitdistrplus::mledist(T_matrix_sum_mid, "lnorm")$estimate
    Lnorm_pval <- 1 - plnorm(TSMatrix, as.numeric(LnormPar[1]), as.numeric(LnormPar[2]))
    return(Lnorm_pval)
  }
  
  # normalize expression matrix
  if(Exp_Norm){
    message("Normalizing the expression matrix")
    ExpMat_SpNorm <- SpMinMax(ExpMat_Sp)
  }
  else{
    ExpMat_SpNorm <- ExpMat_Sp
  }
  Var_wt <- sparseMatrixStats::rowVars(ExpMat_Sp)
  # normalize coordinates
  scalFactor <- 1
  Dimension <- ncol(Coords)
  Coords_Norm_Method <- match.arg(Coords_Norm_Method)
  if(Coords_Norm_Method != "None"){
    message("Normalizing the coords matrix")
    if(Coords_Norm_Method == "Sliced"){
      N_Slice <- 1
      if(Dimension == 3){
        N_Slice <- length(unique(Coords[,3]))
      }
      scalFactor <- exp(mean(log((apply(Coords, 2, quantile, 0.975) - apply(Coords, 2, quantile, 0.025))/0.95)[1:2])) / (nrow(Coords) / N_Slice) ^ (1 / 2)
    }
    if(Coords_Norm_Method == "Overall"){
      scalFactor <- exp(mean(log((apply(Coords, 2, quantile, 0.975) - apply(Coords, 2, quantile, 0.025))/0.95))) / nrow(Coords) ^ (1 / Dimension)
    }
  }
  D1 <- D_1 * scalFactor
  D2 <- D_2 * scalFactor
  
  # calculate variance of local means
  message("Calculating p-values")
  Var_X <- sapply(c(D1, D2), function(DK) VarLocalMeans(DK))
  T_matrix <- Var_X[,2] / Var_X[,1] * Var_wt
  BSP_pval <- Lnorm_estimation(T_matrix)
  
  BSP_outputs <- data.frame(GeneNames = rownames(ExpMat_Sp),
                            P_values = BSP_pval)
  return(BSP_outputs)
}

