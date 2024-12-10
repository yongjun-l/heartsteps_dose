#' Title
#'
#' @param dat data
#' @param sub.id id
#' @param B boot
#'
#' @return ids
#' @importFrom stringr word
#' @export
boot_samples <- function (dat, sub.id, B) 
{
  sub.id <- dat[, sub.id]
  unique.id <- unique(sub.id)
  this.index <- sapply(1:B, function(e) {
    sample(unique.id, size = length(unique.id), replace = TRUE)
  })
  no.repeat.id <- apply(this.index, 2, function(x) {
    temp <- table(x)
    for (i in seq_along(temp)) {
      if (temp[i] > 1) {
        num.appearance <- temp[i]
        x[which(x == names(temp)[i])] <- paste0(names(temp)[i], 
                                                "__", 1:num.appearance)
      }
    }
    x
  })
  output <- apply(no.repeat.id, 2, function(x) {
    temp <- lapply(x, function(x_i) {
      index <- which(sub.id == stringr::word(x_i, 1, sep = "\\__"))
      cbind(index, no.repeat.id = rep(x_i, length(index)))
    })
    dat_return <- Reduce("rbind", temp)
    dat_return <- data.table::as.data.table(dat_return)
    dat_return$index<-as.numeric(dat_return$index)
    # dat_return[, `:=`(index, as.numeric(index))] #Original code
    dat_return
  })
  output
}
