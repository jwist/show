
d <-
  read.csv2("~/Downloads/pvalues_Ranks_controlvseachseverityclass_220222.csv",
            sep = ",")
d <-
  read.csv2("~/Downloads/pvalues_Ranks_controlvseachseverityclass_210922.csv",
            sep = ",")

idx <- match(d$Metabolite.control.vs.C, d$Metabolite.control.vs.B)
d$BC <-
  as.character(factor(
    sign(idx - d$Position.at.control.vs.B),
    labels = c("red", "orange", "green")
  ))
d$idxBC <- idx - d$Position.at.control.vs.B

idx <- match(d$Metabolite.control.vs.D, d$Metabolite.control.vs.C)
d$CD <-
  as.character(factor(
    sign(idx - d$Position.at.control.vs.C),
    labels = c("red", "orange", "green")
  ))
d$idxCD <- idx - d$Position.at.control.vs.C

idx <- match(d$Metabolite.control.vs.E, d$Metabolite.control.vs.D)
d$DE <-
  as.character(factor(
    sign(idx - d$Position.at.control.vs.D),
    labels = c("red", "orange", "green")
  ))
d$idxDE <- idx - d$Position.at.control.vs.D

N <- 1034
d <- d[1:N,]
r <- c(rep("gray75", N), unname(unlist(d[, c("BC", "CD", "DE")])))
sig <- as.numeric(unname(unlist(d[, c(3, 6, 9, 12)])))
sig <- sig < 0.05 # make sig boolean

# contatenating all p-values
p <- unname(unlist(d[, c(2, 5, 8, 11) + 1]))

# concatenating all names
d <- unname(unlist(d[, c(2, 5, 8, 11)]))

n <- t <- f <- d

# replacing names with colors according to class
lipo <- getLipoTable()
idx <- which(d %in% lipo$ID)
d[idx] <- "#0c7209"#"#8888ff"
t[idx] <- "x"
f[idx] <- "2"
all <- idx

# lmwm <- unname(unlist(read.table("sm.name", sep = ",")))
lmwm <- c('tryptophan','3-hydroxykynurenine','3-hydroxyanthranilic acid','kynurenic acid','quinolinic acid','picolinic acid','xanthurenic acid','kynurenine','indole-3-acetic acid','5-hydroxyindole acetic acid','neopterin','serotonin','1-methylhistidine','3-methylhistidine','alanine','alpha-aminobutyric acid','arginine','asparagine','aspartic acid','citrulline','glutamic acid','glutamine','glycine','histidine','isoleucine','leucine','lysine','methionine','ornithine','phenylalanine','proline','serine','taurine','threonine','tyrosine','valine',
          'alanine','creatine','creatinine','glycine','histidine','isoleucine','leucine','lysine','methionine','phenylalanine','tyrosine','valine','acetic acid','citric acid','formic acid','lactic acid','3 hydroxybutyric acid','acetoacetic acid','acetone','pyruvic acid','glucose',
          'serotonin.tryptophan', 'X3.hydroxykynurenine.tryptophan',
          'X3.methylhistidine', 'acetic.acid','Acetoacetic.acid','alpha.aminobutyric.acid','asp...glu...asn...gln','Aspartic.acid','Branched.AA...taurine','Citric.acid','Fisher.s.ratio','Formic.acid','Glutamic.acid','indole.3.acetic.acid','kynurenic.acid','kynurenic.acid.tryptophan','kynurenine.tryptophan','lac.pyr','Lactic.acid','neopterin.tryptophan','picolinic.acid','Pyruvic.acid','quinolinic.acid','quinolinic.acid.tryptophan','serotonin.tryptophan','X1.methylhistidine','X3.hydroxyanthranilic.acid','X3.Hydroxybutyric.acid','X3.hydroxykynurenine','X3.hydroxykynurenine.tryptophan','X3.methylhistidine','X5.hydroxyindole.acetic.acid','xanthurenic.acid',
          'formic.acid', 'pyruvic.acid')
idx <- which(tolower(d) %in% tolower(lmwm))
d[idx] <- "#a31a96"
t[idx] <- "y"
f[idx] <- "3"
all <- c(all, idx)


# lipidClass <- c('CE','CER','DAG','DCER','FFA','HCER','LCER',
#                 'LPC','LPE','LPG','LPI','MAG','PC','PE','PG',
#                 'PI','PS','SM','TAG')
# idx <- which(substr(d, 1, 2) %in% substr(lipidClass, 1, 2))
# d[idx] <- "#000000"
# t[idx] <- sapply(strsplit(t[idx], "\\."), "[", 1)
# f[idx] <- "1"
d[-all] <- "#000000"
t[-all] <- sapply(strsplit(t[-all], "\\."), "[", 1)
f[-all] <- "1"

# create df and sort
df <-
  data.frame(d, r, sig, t, p, n, f, g = rep(c("AB", "CC", "CD", "CE"), each = N))
df %>% group_by(g) %>%
  arrange(t, n, .by_group = TRUE) -> df

# set non significant variable to white and sort by p-value
df$d[!df$sig] <- "#ffffff"
# sorting colors for percent bar
df$s <- (df %>% group_by(g) %>% arrange(f, desc(sig), .by_group = TRUE) %>% select(d))$d

#print summary
df %>% group_by(g, f) %>% summarise(sum(sig)/n())
df %>% group_by(f) %>% summarise(sum(sig))
df %>% filter(sig) %>% arrange(f) -> df2
intersect(df2$n[df2$g == "AB"], df2$n[df2$g == "CC"])

ndx <- setdiff(df2$n[df2$g == "AB"], df2$n[df2$g == "CC"])
tab <- df2[match(ndx, df2$n),]
write.table(tab, file = "notInCCbutinAB.csv", sep = ",")

ndx <- setdiff(df2$n[df2$g == "CC"], df2$n[df2$g == "AB"])
tab <- df2[match(ndx, df2$n),]
write.table(tab, file = "notInABbutinCC.csv", sep = ",")

ndx <- setdiff(df2$n[df2$g == "CC"], df2$n[df2$g == "CD"])
tab <- df2[match(ndx, df2$n),]
write.table(tab, file = "notInCDbutinCC.csv", sep = ",")

ndx <- setdiff(df2$n[df2$g == "CD"], df2$n[df2$g == "CC"])
tab <- df2[match(ndx, df2$n),]
write.table(tab, file = "notInCCbutinCD.csv", sep = ",")

ndx <- setdiff(df2$n[df2$g == "CD"], df2$n[df2$g == "CE"])
tab <- df2[match(ndx, df2$n),]
write.table(tab, file = "notInCEbutinCD.csv", sep = ",")

ndx <- setdiff(df2$n[df2$g == "CE"], df2$n[df2$g == "CD"])
tab <- df2[match(ndx, df2$n),]
write.table(tab, file = "notInCDbutinCE.csv", sep = ",")

# df$s[!df$sig] <- "zz"
# df$s <- c(sort(df$s[1:N]),
#           sort(df$s[(N + 1):(N * 2)]),
#           sort(df$s[(N * 2 + 1):(N * 3)]),
#           sort(df$s[(N * 3 + 1):(N * 4)]))
# df$s[df$s == "zz"] <- "#ffffff"

df$f[df$f == "3"] <- "#a31a96"
df$f[df$f == "2"] <- "#0c7209"
df$f[df$f == "1"] <- "black"

# df$t[df$t == "x"] <- "Lipo"
# df$t[df$t == "y"] <- "LMWM"
svglite::svglite(filename = "severity2.svg", width = 15, height = 4, bg = "transparent")
plotColoredMatrix(
  df$d,
  options = list(
    r = df$f,
    s = df$s,
    sig = df$sig,
    ncol = N,
    by = "row",
    title = "Variables ordered by analytical techniques",
    # title = "Variables ranked by significance (descending from left to right)",
    headers = c(
      "control\nvs.\ngroup B",
      "control\nvs.\ngroup C",
      "control\nvs.\ngroup D",
      "control\nvs.\ngroup E"
    ),
    legend = list(
      color = c("#0c7209",
                "#a31a96",
                "#000000"),
      text = c(
        "IVDr LMWM + AA and Try t-MS assays",
        "IVDr Lipoproteins",
        "MS Lipids"
      )
    )
  )
)
dev.off()



#' print colored matrix
#' @param vec - a vector with colors
#' @param options - options
#' @importFrom grid viewport pushViewport grid.rect grid.newpage
plotColoredMatrix <- function(vec, options = list()) {
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
      stop("fusion >> plotColoredMatrix >> wrong value for 'by' option")
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
  names(ta)[20:21] <- c("LIPO", "LMWM")
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

# plotColoredMatrix(rep("gray75", 11), options = list(ncol = 3, by = "col"))
# plotColoredMatrix(c("#444444", "#444444", "#666666", "#666666"), options = list(ncol = 2, by = "row"))
# plotColoredMatrix(c("#444444", "#444444", "#666666", "#666666"), options = list(ncol = 2, by = "col"))
# plotColoredMatrix(rep("gray75", 11), options = list(by = "row"))
