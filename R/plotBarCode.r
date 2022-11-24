library(tidyverse)
library(grid)
d <- read.csv2("~/Downloads/lipid_paper_barcode_plot_data.csv",
            sep = ",")
d$class <- getLipidsInfo(d$lipid)$nC
d %>% gather(key, value, impute:concentration_signal_correct, -lipid) -> df

df$sig <- !is.na(df$value)

lipidClass <- getLipidsInfo(df$lipid)

df$class <- lipidClass$class
classes <- unique(df$class)
classColors <- viridis::viridis(length(classes))
df$d <- "black"
for (i in 1:length(classes)) {
  df$d[df$class == classes[i]] <- classColors[i]
}
df$t <- df$class


N <- nrow(d)
df$g = rep(c("AB", "CC", "CD"), each = N)

df <- df %>% group_by(g) %>% arrange(class, .by_group = TRUE)

# set non significant variable to white and sort by p-value
df$d[!df$sig] <- "#ffffff"
# sorting colors for percent bar
df$s <- (df %>% group_by(g) %>% arrange(class, desc(sig), .by_group = TRUE) %>% select(g,d))$d
df$r <- (df %>% group_by(g) %>% arrange( desc(sig), class, .by_group = TRUE) %>% select(g,d))$d

svglite::svglite(filename = "lipidsBarcodeByClass.svg", width = 25, height = 4, bg = "transparent")
plotBarcode(
  df$d,
  options = list(
    r = df$r,
    s = df$s,
    sig = !is.na(df$value),
    ncol = nrow(d),
    by = "row",
    title = "Variables ordered by analytical techniques",
    # title = "Variables ranked by significance (descending from left to right)",
    headers = c(
      "impute",
      "concentration",
      "concentration with signal correction"
    ),
    legend = list(
      color = classColors,
      text = classes
    )
  )
)
dev.off()



#' print colored matrix
#' @param vec - a vector with colors
#' @param options - options
#' @importFrom grid viewport pushViewport grid.rect grid.newpage
plotBarcode <- function(vec, options = list()) {
  if ("ncol" %in% names(options)) {
    ncol <- options$ncol
  } else {
    ncol <- 1
  }

  if ("by" %in% names(options)) {
    by <- options$by
  } else {
    by <- "col"
  }

  if ("headers" %in% names(options)) {
    headers <- options$headers
  } else {
    headers <- rep("", ncol)
  }

  if ("legend" %in% names(options)) {
    legend <- options$legend
  } else {
    legend <- rep("", length(unique(vec)))
  }

  if ("r" %in% names(options)) {
    r <- options$r
  } else {
    r <- rep("gray75", length(vec))
  }

  if ("s" %in% names(options)) {
    s <- options$s
  } else {
    s <- rep("gray75", length(vec))
  }

  if ("sig" %in% names(options)) {
    sig <- options$sig
  } else {
    sig <- rep(TRUE, length(vec))
  }

  if ("title" %in% names(options)) {
    title <- options$title
  } else {
    title <- ""
  }


  l <- length(vec)
  # create
  grid.newpage()
  makeGrid(l, ncol)
  # grid.rect()

  # fill colors
  h <- ceiling(l / ncol)
  for (i in 1:l) {
    if (by == "col") {
      rdx <- rep(c(0:(h - 1)), ncol)
      ndx <- rep(c(0:(ncol - 1)), each = h)
    } else if (by == "row") {
      rdx <- rep(c(0:(h - 1)), each = ncol)
      ndx <- rep(c(0:(ncol - 1)), h)
    } else {
      stop("show >> plotBarcode >> wrong value for 'by' option")
    }

    vp <- viewport(
      x = unit(0.5 + ndx[i], "native"),
      y = unit(h - 0.5 - rdx[i], "native"),
      width = unit(1, "native"),
      height = unit(0.9, "native")
    )
    pushViewport(vp)
    if (sig[i]) {
      grid.rect(
        x = unit(0.5, "native"),
        y = unit(0.6, "native"),
        height = unit(0.6, "native"),
        gp = gpar(fill = vec[i], col = NA)
      )

    }
    grid.rect(
      x = unit(0.5, "native"),
      y = unit(0.95, "native"),
      height = unit(0.05, "native"),
      gp = gpar(fill = r[i], col = NA)
    )
    grid.rect(
      x = unit(0.5, "native"),
      y = unit(0.1, "native"),
      height = unit(0.25, "native"),
      gp = gpar(fill = s[i], col = NA)
    )
    upViewport()
  }
  upViewport()

  # print headers
  if (FALSE) {
    vp <- viewport(
      x = unit(0.4, "npc"),
      y = unit(0.975, "npc"),
      width = unit(0.6, "npc"),
      height = unit(0.05, "npc"),
      xscale = c(0, ncol),
      yscale = c(0, 1)
    )
    pushViewport(vp)
    # grid.rect()
    for (i in 1:ncol) {
      vp <- viewport(
        x = unit(i - 1 + 0.5, "native"),
        y = unit(0.5, "native"),
        width = unit(0.8, "native"),
        height = unit(0.8, "native")
      )
      pushViewport(vp)
      # grid.rect()
      grid.text(
        headers[i],
        gp = gpar(cex = 0.7, alpha = 0.6),
        just = "center",
        x = unit(0.5, "native"),
        y = unit(0.5, "native")
      )
      upViewport()
    }
    upViewport()
  } else {
    vp <- viewport(
      x = unit(0.05, "npc"),
      y = unit(0.5, "npc"),
      width = unit(0.1, "npc"),
      height = unit(0.9, "npc"),
      xscale = c(0, 1),
      yscale = c(0, h)
    )
    pushViewport(vp)
    # grid.rect()
    for (i in 1:h) {
      vp <- viewport(
        x = unit(0.5, "native"),
        y = unit(h + 0.5 - i, "native"),
        width = unit(0.8, "native"),
        height = unit(0.8, "native")
      )
      pushViewport(vp)
      # grid.rect()
      grid.text(
        headers[i],
        gp = gpar(cex = 0.7, alpha = 0.6),
        just = "center",
        x = unit(0.5, "native"),
        y = unit(0.5, "native")
      )
      upViewport()
    }
    upViewport()
  }

  # print legend
  ll <- length(legend$color)
  vp <- viewport(
    x = unit(0.85, "npc"),
    y = unit(0.5, "npc"),
    width = unit(0.3, "npc"),
    height = unit(0.2, "npc"),
    xscale = c(0, 2),
    yscale = c(0, ll)
  )
  pushViewport(vp)
  # grid.rect()
  for (i in 1:ll) {
    vp <- viewport(
      x = unit(0.5, "native"),
      y = unit(ll + 0.5 - i, "native"),
      width = unit(0.8, "native"),
      height = unit(1, "native")
    )
    pushViewport(vp)
    grid.rect(
      x = unit(0.15, "native"),
      width = unit(0.3, "native"),
      height = unit(0.8, "native"),
      gp = gpar(fill = legend$color[i], col = NA)
    )
    grid.text(
      legend$text[i],
      gp = gpar(cex = 0.7, alpha = 0.6),
      just = "left",
      x = unit(0.4, "native"),
      y = unit(0.5, "native")
    )
    upViewport()
  }
  upViewport()

  # print xlab
  vp <- viewport(
    x = unit(0.4, "npc"),
    y = unit(0.025, "npc"),
    width = unit(0.6, "npc"),
    height = unit(0.05, "npc"),
    xscale = c(0, 1),
    yscale = c(0, 1)
  )
  pushViewport(vp)
  grid.text(
    title,
    gp = gpar(cex = 0.7, alpha = 0.6),
    just = "center",
    x = unit(0.5, "native"),
    y = unit(0.7, "native")
  )
  upViewport()

  # labels
  vp <- viewport(
    x = unit(0.4, "npc"),
    y = unit(0.94, "npc"),
    width = unit(0.6, "npc"),
    height = unit(0.05, "npc"),
    xscale = c(0, ncol),
    yscale = c(0, 1)
  )
  pushViewport(vp)
  # grid.rect()
  ta <- table(df$t[1:ncol])
  cu <- cumsum(ta)
  ce <- cu - ta / 2
  for (i in 1:length(ta)) {
    vp <- viewport(
      x = unit(ce[i] + 0.5, "native"),
      y = unit(0.5, "native"),
      width = unit(ta[i], "native"),
      height = unit(0.8, "native")
    )
    pushViewport(vp)
    # grid.rect()
    grid.lines(x = unit(c(0, 1), "native"),
               y = unit(c(0, 0), "native"))
    grid.lines(
      x = unit(c(0, 0), "native"),
      y = unit(c(-21, 0.2), "native"),
      gp = gpar(lty = 3, alpha = 0.3)
    )
    grid.lines(
      x = unit(c(1, 1), "native"),
      y = unit(c(-21, 0.2), "native"),
      gp = gpar(lty = 3, alpha = 0.3)
    )
    grid.text(
      names(ta)[i],
      gp = gpar(cex = 0.5, alpha = 0.6),
      just = "left",
      x = unit(0.5, "native"),
      y = unit(0.5, "native"),
      rot = 90
    )
    upViewport()
  }
  upViewport()
}

#' create grid with grid package
#' @param length - the title
#' @param ncol - position
#' @importFrom grid viewport pushViewport
makeGrid <- function(length, ncol) {
  vp <- viewport(
    x = unit(0.4, "npc"),
    #normalized unit (0 -> 1)
    width = 0.6,
    height = 0.85,
    xscale = c(0, ncol),
    yscale = c(0, ceiling(length / ncol))
  )
  pushViewport(vp)
}

getLipidsInfo <- function(lipid) {
  lmc <- strsplit(lipid, "\\(")
  lmc <- rapply(lmc, function(x) c(x[1], gsub("\\)", "", x[2])), how = "replace")
  r <- list()
  for (i in 1:length(lmc)) {
    struc <- strsplit(lmc[[i]][2], "_")[[1]]
    if (length(struc) == 1 | lmc[[i]][1] == "TAG") {
      totalCarbon <- gsub("[a-zA-Z\\-]+",
                          "",
                          strsplit(struc, ":")[[1]][1])
      if (lmc[[i]][1] == "TAG") {
        sideChain <- gsub("[a-zA-Z\\-]+",
                          "",
                          strsplit(lmc[[i]][2], "_")[[1]][2])
      } else {
        unsat <- gsub("[a-zA-Z\\-]+",
                      "",
                      strsplit(struc, ":")[[1]][2])
        sideChain <- paste0(totalCarbon, ":", unsat)
      }
    } else {
      t <- strsplit(struc, ":")
      sc <- unlist(lapply(t, function(x)
        as.numeric(
          gsub("[a-zA-Z\\-]+",
               "",
               x[1])
        )))
      unsatSc <- unlist(lapply(t, function(x)
        as.numeric(
          gsub("[a-zA-Z\\-]+",
               "",
               x[2])
        )))
      totalCarbon <- sum(sc)
      sideChain <- paste0(sc, ":", unsatSc)
      unsat <- sum(unsatSc)
    }
    r[[i]] <- c(lmc[[i]][1], totalCarbon, unsat, sideChain)
  }
  l <- max(unlist(lapply(r, function(x) length(x))))
  r <- lapply(r, function(x) c(x, rep(NA, l-length(x))))
  lipidClass <- data.frame(do.call("rbind", r))
  colnames(lipidClass) <- c("class", "nC", "r", "sc1", "sc2")
  return(lipidClass)
}
# plotColoredMatrix(rep("gray75", 11), options = list(ncol = 3, by = "col"))
# plotColoredMatrix(c("#444444", "#444444", "#666666", "#666666"), options = list(ncol = 2, by = "row"))
# plotColoredMatrix(c("#444444", "#444444", "#666666", "#666666"), options = list(ncol = 2, by = "col"))
# plotColoredMatrix(rep("gray75", 11), options = list(by = "row"))
