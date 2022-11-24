
library(networkD3)

# Load energy projection data
URL <- "https://cdn.rawgit.com/christophergandrud/networkD3/master/JSONdata/energy.json"
Energy <- jsonlite::fromJSON(URL)


# Now we have 2 data frames: a 'links' data frame with 3 columns (from, to, value), and a 'nodes' data frame that gives the name of each node.
head( Energy$links )
head( Energy$nodes )
# Thus we can plot it
p <- sankeyNetwork(Links = Energy$links, Nodes = Energy$nodes, Source = "source",
                   Target = "target", Value = "value", NodeID = "name",
                   units = "f", fontSize = 15, nodeWidth = 100)
p


df <- data.frame(
  links = I(data.frame(
    source = c(0,1,2,3,4,5,5),
    target = rep(6,7),
    value = c(2095,5029,4248,1500, 19114, 900, 1))),
  nodes = I(data.frame(name = c("RAINE", "BHAS", "HIMS", "PLSAW", "ASPREE", "SCOT HEART", ""))))
p <- sankeyNetwork(Links = df$links, Nodes = df$nodes, Source = "source",
              Target = "target", Value = "value", NodeID = "name",
              units = "", fontSize = 15, nodeWidth = 10, height = 400, width = 150);
svglite::svglite("cohortIn.svg", bg = "transparent", width = 5, height = 10)
print(p)
dev.off()


nodes = data.frame("name" =
                     c("Node A",
                       "Node B",
                       "Node C",
                       "Node D"))
links = as.data.frame(matrix(c(
  0, 1, 10, # Each row represents a link. The first number
  1, 2, 20, # represents the node being conntected from.
  1, 3, 30, # the second number represents the node connected to.
  2, 3, 40),# The third number is the value of the node
  byrow = TRUE, ncol = 3))
names(links) = c("source", "target", "value")
sankeyNetwork(Links = links, Nodes = nodes,
              Source = "source", Target = "target",
              Value = "value", NodeID = "name",
              fontSize= 12, nodeWidth = 30)
