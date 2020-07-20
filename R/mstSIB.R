#' The mstSIB test for MSTs
#'
#' This function allows the detection of itemwise DIF using the mstSIB test.
#'
#' Author: Mark J. Gierl, with minor changes by Rudolf Debelak and Dries Debeer
#'
#' @param resp A data frame containing the response matrix. Rows correspond to respondents, columns to items.
#' @param DIF_covariate A vector indicating the membership to the reference (0) and focal (1) groups.
#' @param theta A vector of ability estimates for each respondent.
#' @param see A vector of the standard error of the ability estimates for each respondent.
#' @param NCell The initial number of cells for estimating the overall ability difference between the focal and reference groups.
#' @param cellmin Minimum number of respondents per cell for the focal and reference group. Cells with fewer respondents are discarded.
#' @param pctmin Minimum rate of focal and reference group that should be used for estimating the over ability difference between focal and groups after discarding cells with few respondents.
#'
#' @return A list with four elements. The first element is the response matrix, the second element is the name of
#' the DIF covariate, and the third element is the name of the test. The fourth element is a matrix where each
#' row corresponds to an item. The columns correspond to the following entries:
#' \describe{
#'   \item{Beta}{The estimated weighted ability difference between the focal and reference groups.}
#'   \item{Vars}{The estimation error of the weighted ability difference between the focal and reference groups.}
#'   \item{N_R}{The number of respondents in the reference group.}
#'   \item{N_F}{The number of respondents in the focal group.}
#'   \item{NCell}{The initial number of cells for estimating the overall ability
#'   difference between the focal and reference groups.}
#'   \item{p_value}{The p-value of the null hypothesis that the ability difference
#'   between the focal and reference groups is 0.}
#'   }
#'
#' @examples 
#' data("toydata")
#' resp <- toydata$resp
#' group_categ <- toydata$group_categ
#' theta_est <- toydata$theta_est
#' see_est <- toydata$see_est
#' mstSIB(resp = as.data.frame(resp), theta = theta_est,
#' DIF_covariate = group_categ, see = see_est)
#'
#' @export

mstSIB <- function(resp, DIF_covariate, theta = NULL, see = NULL,
                   cellmin = 3, pctmin = .9, NCell = 80){

  # get call
  call <- match.call()

  # a theta-argument is required
  if(is.null(theta)) stop("'theta'-argument is missing. Include a vector with the estimated 'theta'-values.", call. = FALSE)

  # a see-argument is required
  if(is.null(see)) stop("'see'-argument is missing. Include a vector with the estimated standard errors of the 'theta'-values.", call. = FALSE)


  # get the DIF_covariate name
  DIF_covariate_name <- as.character(deparse(call$DIF_covariate))

  # only works for two groups coded 0 and 1
  DIF_covariate <- as.factor(DIF_covariate)
  stopifnot(nlevels(DIF_covariate) == 2)
  levels(DIF_covariate) <- c(0, 1)
  DIF_covariate <- as.numeric(as.character(DIF_covariate))


  # get number of items
  nItem <- ncol(resp)

  # get/set item names
  colnames(resp) <- itemnames <- 'if'(is.null(colnames(resp)),
                                      sprintf(paste("it%0", nchar(nItem),
                                                    "d", sep=''),
                                              seq_len(nItem)),
                                      colnames(resp))

  ##insert by variable
  Sif <- cbind(theta, see, DIF_covariate, resp)
  BetaOut<-matrix(numeric(0),dim(Sif)[2]-3,6)
  rownames(BetaOut) <- colnames(resp)
  colnames(BetaOut) <- c("stat", "SE", "N_R", "N_F", "NCell", "p_value")
  ##Start here
  for(inum in 4:dim(Sif)[2]){
    Rif <- Sif[Sif[, 3] == 0 & !is.na(Sif[, inum]), ]
    Fif <- Sif[Sif[, 3] == 1 & !is.na(Sif[, inum]), ]
    if(nrow(Rif) > 0 & nrow(Fif) > 0){
      RSR <- Rif[, inum]
      FSR <- Fif[, inum]

      ## Splitting the file into the focus and reference group, their item response starts at col 6
      EThetaF <- Fif[1:dim(Fif)[1],1]
      ESEF<-Fif[1:dim(Fif)[1],2]
      EThetaR<-Rif[1:dim(Rif)[1],1]
      ESER<-Rif[1:dim(Rif)[1],2]
      MeanR<-mean(EThetaR)
      VarR<-var(EThetaR)
      MeanF<-mean(EThetaF)
      VarF<-var(EThetaF)
      TMin<-min(min(EThetaR),min(EThetaF))
      TMax<-max(max(EThetaR),max(EThetaF))
      RI<-(ESER^-2)
      FI<-(ESEF^-2)
      MeanRI<-mean(RI)
      MeanFI<-mean(FI)

      ## fixed, the equation is be = rmean + (1. - (1./rinfomean)/rvar)*(thetaro(j) - rmean)
      ## we need test information
      AThetaR<-MeanR+(1-(1/MeanRI/VarR))*(EThetaR - MeanR)
      AThetaF<-MeanF+(1-(1/MeanFI/VarF))*(EThetaF - MeanF)

      ## Sort all examinees and their responses by their thetahats
      RefMain<-cbind(AThetaR,RSR)
      FocMain<-cbind(AThetaF,FSR)

      RefMain<-RefMain[order(RefMain[,1]),]
      FocMain<-FocMain[order(FocMain[,1]),]

      ## Defining Min and Max for interval establishment
      TMin<-min(min(RefMain[,1]),min(FocMain[,1]))
      TMax<-max(max(RefMain[,1]),max(FocMain[,1]))

      ##Define Initial number of cells
      ##try here first, finding and counting for bins
      RefInt<-findInterval(RefMain[,1],(TMin+((TMax-TMin)/NCell)*0:NCell), rightmost.closed = FALSE, all.inside = FALSE)
      FocInt<-findInterval(FocMain[,1],(TMin+((TMax-TMin)/NCell)*0:NCell), rightmost.closed = FALSE, all.inside = FALSE)
      CellCountR<-0
      CellCountF<-0
      for(i in 1:(NCell+1)){
        CellCountR[i]<-length(RefInt[RefInt==i])
        CellCountF[i]<-length(FocInt[FocInt==i])
        if ((CellCountR[i]<cellmin)||(CellCountF[i]<cellmin)){
          CellCountR[i]<-0
          CellCountF[i]<-0
          RefInt[RefInt==i]<-0
          FocInt[FocInt==i]<-0
        }
      }
      while(((sum(CellCountR)<=pctmin*dim(RefMain)[1])&&(NCell>5))||((sum(CellCountF)<=pctmin*dim(FocMain)[1])&&(NCell>5))||((sum(CellCountR)+sum(CellCountF))<=(pctmin*dim(RefMain)[1]+pctmin*dim(FocMain)[1])&&(NCell>5))) {
        NCell<-NCell-4
        RefInt<-findInterval(RefMain[,1],(TMin+((TMax-TMin)/NCell)*0:NCell), rightmost.closed = FALSE, all.inside = FALSE)
        FocInt<-findInterval(FocMain[,1],(TMin+((TMax-TMin)/NCell)*0:NCell), rightmost.closed = FALSE, all.inside = FALSE)
        CellCountR<-0
        CellCountF<-0
        for(i in 1:(NCell+1)){
          CellCountR[i]<-length(RefInt[RefInt==i])
          CellCountF[i]<-length(FocInt[FocInt==i])
          if ((CellCountR[i]<cellmin)||(CellCountF[i]<cellmin)){
            CellCountR[i]<-0
            CellCountF[i]<-0
            RefInt[RefInt==i]<-0
            FocInt[FocInt==i]<-0
          }
        }
      }
      ##check numbers for bins
      ##NCell
      ##sum(CellCountF)
      ##sum(CellCountR>1)
      ##CellCountF[1]
      ##CellCountR[1]
      ##RefInt[RefInt==1]
      ##plot(RefInt)

      ##item proportion for bins
      beta<-0
      items<- if(is.null(dim(RSR))) 1 else (dim(RSR)[2])
      vars<-0
      for(j in 1:items){
        ybarR<-0
        ybarF<-0
        uf2sum<-0
        #ufsum<-0
        #ursum<-0
        ur2sum<-0
        for(i in 1:NCell+1){
          uf2sum[i]<-sum(FocMain[FocInt==i,j+1])^2
          #ufsum[i]<-sum(FocMain[FocInt==i,j+1])
          #ursum[i]<-sum(RefMain[RefInt==i,j+1])
          ur2sum[i]<-sum(RefMain[RefInt==i,j+1])^2
          ybarR[i]<-sum(RefMain[RefInt==i,j+1])/CellCountR[i]
          ybarF[i]<-sum(FocMain[FocInt==i,j+1])/CellCountF[i]
        }
        wt<-(CellCountR+CellCountF)/(sum(CellCountR)+sum(CellCountF))
        wtsum<-sum(wt)
        varr<-(ur2sum-(CellCountR*ybarR*ybarR))/(CellCountR-1)
        varf<-(uf2sum-(CellCountF*ybarF*ybarF))/(CellCountF-1)
        varr
        varf
        bbg<-(ybarR-ybarF)*wt
        var<-wt*wt*((1/CellCountR)*varr + (1/CellCountF)*varf)
        ##plot((PropCorR/CellCountR)-(PropCorF/CellCountF))
        beta[j]<-sum(bbg, na.rm=TRUE)
        vars[j]<-sum(var,na.rm=TRUE)
      }
    } else{beta<- NA
    vars <- NA
    NCell <- NA}
    BetaOut[inum-3,1]<-beta
    BetaOut[inum-3,2]<-vars
    BetaOut[inum-3,3]<-dim(Rif)[1]
    BetaOut[inum-3,4]<-dim(Fif)[1]
    BetaOut[inum-3,5]<-NCell
    BetaOut[inum-3,6]<-2 * stats::pnorm(-abs(beta/vars))

  }
  return(list(resp = resp,
              DIF_covariate = DIF_covariate_name,
              test = "SIB-test",
              results = as.data.frame(BetaOut)))
}

