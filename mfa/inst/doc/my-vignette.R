## ------------------------------------------------------------------------
library(mfa)
df <- read.csv("../data/wines.csv", stringsAsFactors = FALSE)
sets=list(2:7,8:13,14:19)
mfa1 <- mfa(df, sets, ncomps = NULL, center = TRUE, scale = TRUE)

## ------------------------------------------------------------------------
plot_comp(mfa1)

## ------------------------------------------------------------------------
plot_pfs(mfa1)

## ------------------------------------------------------------------------
graphics.off()
plot_loadings(mfa1)

## ------------------------------------------------------------------------
eigentbl(mfa1)

## ------------------------------------------------------------------------
ob2dim(mfa1)

## ------------------------------------------------------------------------
var2dim(mfa1)

## ------------------------------------------------------------------------
tbl2dim(mfa1)

## ------------------------------------------------------------------------
RV(df[,2:7],df[8:13])

## ------------------------------------------------------------------------
RVtbl(df,sets=list(2:7,8:13,14:19))

