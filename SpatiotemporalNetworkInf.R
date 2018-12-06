setwd("C:/Users/fbx5002/Desktop/ThesisCodes")

# Solution for the permit of Rtools
Sys.setenv(PATH = paste(Sys.getenv("PATH"), "D:/Rtools/bin/",
                        "D:/Rtools/mingw_64/bin", sep = ";")) #for 64 bit version
Sys.setenv(BINPREF = "D:/Rtools/mingw_64/bin")
assignInNamespace("version_info", c(devtools:::version_info, list("3.5" = list(version_min = "3.3.0", version_max = "99.99.99", path = "bin"))), "devtools")
find_rtools()
# print(tools::showNonASCIIfile("C:/Users/fbx5002/Desktop/ThesisCodes/thesis.R"))
# rsconnect::setAccountInfo(name='fangcaoxu',
#                           token='DC2393B3873096BE30B023787D8B51C5',
#                           secret='OeeTaz5KGj7cQ9hKFi9qVXXxUD0i1awzQLO71CAQ')

# Load Libraries ---------------------------------------------------------------------------------------
#library(BayesX)
library(combinat)
library(coda)
library(dplyr)
library(fields)
library(ggplot2)
library(gridExtra)
#library(grid)
library(igraph)
library(leaflet)
library(maptools)
library(mapview)
#library(networkSpatial)
#library(NetworkInference, lib.loc="D:/R-3.5.1/library")
#library(netdiffuseR)
library(plotly)
library(plyr)
library(raster)
library(rsconnect)
library(reshape2)
library(rgdal)
library(rgeos)
library(RColorBrewer)
#library(R2BayesX)
library(shiny)
library(sp)
#library(spatsurv)
#library(stringi)
library(statnet)
#library(survival)
#library(spBayesSurv)
library(tidyr)
library(tidyverse)

# Functions --------------------------------------------------------------------------------------------
# declared variables in functions are local. If you want to make it global, use <<- or assign() function
### Check Reciprocity
isreciprocal  <- function(d){
  d1.vec <- apply(d, 1, paste, collapse = " ")
  d2.vec <- apply(cbind(d[,2],d[,1]), 1, paste, collapse = " ")
  d1.in.d2.rows <- d[d1.vec %in% d2.vec,]
  if(nrow(d1.in.d2.rows)==0){
    return("There are no two rows being reciprocal with each other!")
  }else{
    d1.in.d2.rows
  }
} 
### Random Network Builder
redgebuilder <- function (nodesid, num_nodes, num_edges){
  edgematrix <- matrix(0, num_nodes,num_nodes)
  lowertri <- which(lower.tri(edgematrix))
  if (num_edges <= length(lowertri)){
    indice <- sample(lowertri, num_edges, replace=FALSE)
    indiceupper <- sample(indice, round(num_edges/2), replace=FALSE)
    indicelower <- indice[!(indice %in% indiceupper)]
    edgematrix[indiceupper]<- 1
    edgematrix <- t(edgematrix)
    edgematrix[indicelower]<- 1
    redgelist <- data.frame(get.edgelist(graph.adjacency(edgematrix)))
    return(redgelist)
  }else{
    return("The edge number has exceeded the maximum allowance for a non-reciprocal directed network.")
  }
}
### Small-world Network Builder
swedgebuilder <- function (nodesid, num_nodes, num_edges){
  k <- round(num_edges/num_nodes)    # neighbors
  swedgelist <- setNames(data.frame(matrix(ncol = 3, nrow = 0)),c("nodes","neinodes","dist"))
  for(i in nodesid){
    if(!(i %in% swedgelist$neinodes)){
      kneighbors<- rownames_to_column(sort(data.frame(dist.mat[i,])),var = "neinodes")[2:(k+1),]
      colnames(kneighbors)[2] <- "dist"
      kneighbors$neinodes<-sapply(kneighbors$neinodes,as.numeric)
      swedgelist<- unique(data.frame(rbind(swedgelist,cbind(nodes=i,kneighbors))))
    }
    else{                           # avoid the network reciprocity
      node <- swedgelist$nodes[which(swedgelist$neinodes==i)]
      kneighbors<- rownames_to_column(sort(data.frame(dist.mat[i,])),var = "neinodes")
      colnames(kneighbors)[2] <- "dist"
      kneighbors$neinodes<-sapply(kneighbors$neinodes,as.numeric)
      kneighbors <- kneighbors[-which(kneighbors$neinodes %in% node),][2:(k+1),]
      swedgelist<- unique(data.frame(sapply(rbind(swedgelist,cbind(nodes=i,kneighbors)),as.numeric)))
    }
  }
  swedgelist
}
### Adjacency Matrix builder
adjacencybuilder <- function(nodesid, edgelist){
  adj.mat <- matrix(0, ncol = length(nodesid), nrow = length(nodesid))
  for(k in 1:nrow(edgelist)) {
    i <- which(nodesid == edgelist[k, 1])
    j <- which(nodesid == edgelist[k, 2])
    adj.mat[i, j] <- 1
  }
  adj.mat
}
### Spatial Lines builder
splinesbuilder <- function(points, edgelist){         # points have three columns: ID, latitude, and longitude
  linelist <- list()
  for(i in 1:nrow(edgelist)){
    line <- Line(rbind(points[edgelist[i,1],c(2,3)], points[edgelist[i,2],c(2,3)]))
    linelist[[i]] <- Lines(line,ID= i)
  }
  spLines <- SpatialLines(linelist)
}
### Simulate Cascades for Given Starting Node over the Network
simulatecascades <- function(start_nodes, g){
  shortestpath <- as.vector(igraph::shortest_paths(g, from = start_nodes, mode = "out", output = "vpath", predecessors=TRUE))
  parentnode <- as.vector(shortestpath$predecessors)
  vpath <- matrix(lapply(shortestpath$vpath,unname))  # if the node name is not same with id here, don't use unname
  vpath[which(is.na(parentnode))] <- NA
  dists <- as.vector(igraph::distances(g, v = start_nodes, mode = "out"))
  tmp <- data.frame("node_name" = V(g)$name, "event_time"= as.numeric(dists),"parent_node" = parentnode, "v_path"=vpath, "cascade_id" = start_nodes,stringsAsFactors=FALSE)
  if(length(which(tmp$event_time==Inf)) != 0){
    tmp <- tmp[-which(tmp$event_time==Inf),]
  }
  if(length(tmp) == 1){                               # if the node is an isolated node
    out <- matrix(vector(), 0, 5) 
  } else{
    tmp <- tmp[order(tmp$event_time),]
    out <- tmp
  }
  out
}
### Transform Long Formatted Cascades to a Formatted List (each row describes one event in the cascade)
as_cascade_long <- function (cascade){
  cascade$node_name <- as.character(cascade$node_name)
  cascade$cascade_id <- as.character(cascade$cascade_id)
  node_names <- unique(cascade$node_name)
  splt <- split(cascade, f = cascade$cascade_id) # split(): return a list
  cascade_nodes <- lapply(splt, function(x) x$node_name)
  cascade_times <- lapply(splt, function(x) as.numeric(x$event_time))
  out <- list("cascade_nodes" = cascade_nodes, "cascade_times" = cascade_times, "node_names" = node_names)
  class(out) <- c("cascade", "list")
  out <- order_cascade(out)
  return(out)
}
### Order Cascade by Time
order_cascade<- function(out) {
  casc_id <- names(out$cascade_times)
  sort_times <- function(x) return(out$cascade_times[[x]][orderings[[x]]])
  sort_nodes <- function(x) return(out$cascade_nodes[[x]][orderings[[x]]])
  orderings <- lapply(out$cascade_times, order)
  times <- lapply(c(1:length(out$cascade_times)), sort_times)
  ids <- lapply(c(1:length(out$cascade_nodes)), sort_nodes)
  out$cascade_nodes <- ids    
  out$cascade_times <- times
  names(out$cascade_nodes) <- casc_id
  names(out$cascade_times) <- casc_id
  return(out)
}
### Plot Spatial Network
plotspatialnetwork <- function (sppoints, splines){
  m <- leaflet() %>% addTiles('http://{s}.basemaps.cartocdn.com/dark_all/{z}/{x}/{y}.png', 
                              attribution='Map tiles by <a href="http://stamen.com">Stamen Design</a>, <a href="http://creativecommons.org/licenses/by/3.0">CC BY 3.0</a> &mdash; 
                              Map data &copy; <a href="http://www.openstreetmap.org/copyright">OpenStreetMap</a>') %>%
    addCircles(data = sppoints, radius=10000, stroke=TRUE, color = "red", opacity = 1, weight = 1, fillOpacity = 1)%>%
    addPolylines(data = splines,color = "white", weight = 0.5) %>%
    setView(lng = -96, lat = 37, zoom =4) 
  m
}

###################################################################################################################################
###########################################      Codes        Start         Here        ###########################################
###########################################      Codes        Start         Here        ###########################################
###########################################      Codes        Start         Here        ###########################################
###################################################################################################################################


# Section 1: Simulate Mutiple Trees --------------------------------------------------
set.seed(111)
samp.nodes <- data.frame(nodes = c("Chicago", "Pittsburgh", "Washington", "Philadelphia", "New York", "Boston"))
samp.nodes$x <- c(-87.6298,-79.9959,-77.0369,-75.1652,-74.0060,-71.0589)
samp.nodes$y <- c(41.8781, 40.4406, 38.9072,39.9526,40.7128,42.3601)
samp.nodes$x <- (samp.nodes$x +88)/8
samp.nodes$y <- (samp.nodes$y -38)/8
samp.edges<- rbind(c("Washington","Philadelphia"),c("Washington","Pittsburgh"),c("Washington","Chicago"),
              c("Philadelphia","Pittsburgh"), c("Pittsburgh","Chicago"), 
              c("Philadelphia", "New York"),c("Boston","New York"),c("Chicago","New York"),
              c("Pittsburgh","New York"),c("Chicago","Boston"))
samp.g <- graph.data.frame(samp.edges,samp.nodes,directed = TRUE)
par(mfrow=c(1,3))
E(samp.g)[c(1,4,5)]$color = "red"
E(samp.g)[c(1,4,5)]$width = 2
E(samp.g)[-c(1,4,5)]$color = "black"
E(samp.g)[-c(1,4,5)]$width = 1
plot(samp.g, layout=as.matrix(samp.nodes[,c(2,3)]), vertex.color= "white",
     vertex.size=30, vertex.label.color="black", vertex.label.cex=1.5, 
     vertex.frame.color ="black", edge.arrow.size=.1)
E(samp.g)[c(1,2,5)]$color = "red"
E(samp.g)[c(1,2,5)]$width = 2
E(samp.g)[-c(1,2,5)]$color = "black"
E(samp.g)[-c(1,2,5)]$width = 1
plot(samp.g, layout=as.matrix(samp.nodes[,c(2,3)]), vertex.color= "white",
     vertex.size=30, vertex.label.color="black", vertex.label.cex=1.5, 
     vertex.frame.color ="black",edge.arrow.size=.1)
E(samp.g)[c(1,3,4)]$color = "red"
E(samp.g)[c(1,3,4)]$width = 2
E(samp.g)[-c(1,3,4)]$color = "black"
E(samp.g)[-c(1,3,4)]$width = 1
plot(samp.g, layout=as.matrix(samp.nodes[,c(2,3)]), vertex.color= "white",
     vertex.size=30, vertex.label.color="black", vertex.label.cex=1.5,
     vertex.frame.color ="black", edge.arrow.size=.1)

# Section 2: Simulate Geographic Network Nodes ---------------------------------------------
num_nodes <- 100
num_edges <- 600
crdref <- CRS('+proj=longlat +datum=WGS84')
extent <- extent(-124.7844079, -66.9513812, 24.7433195, 49.3457868) # extent follows xmin, xmax, ymin, ymax
map <- crop(spTransform(readOGR("cb_2017_us_state_5m","cb_2017_us_state_5m"),crdref),extent)
sppoints<-spsample(map, n=num_nodes, type= "random", proj4string=crdref) # create random spatial points
id <- data.frame(id = 1:num_nodes)
sppoints.df <- SpatialPointsDataFrame(sppoints,id)
points<- data.frame(sppoints.df)
points$optional <- NULL
dist.mat <- spDists(sppoints)     # create distance matrix
dist.mat <- round(dist.mat*1000)  # in meters

# Section 3: Simulate Edges of Random / Small-world Network -----------------------------------------------
redgelist <- redgebuilder(points$id, num_nodes, num_edges)
isreciprocal(redgelist)
swedgelist <- swedgebuilder(points$id, num_nodes, num_edges)
isreciprocal(swedgelist[,c(1,2)])
rspLines <- splinesbuilder(points, redgelist)
rspLines.df <- SpatialLinesDataFrame(rspLines, redgelist)
swspLines <- splinesbuilder(points, swedgelist)
swspLines.df <- SpatialLinesDataFrame(swspLines, swedgelist)
plotspatialnetwork(sppoints, rspLines)
#mapshot(..., file = "RandomNetwork.jpg")
plotspatialnetwork(sppoints, swspLines)
#mapshot(..., file = "SmallWorldNetwork.jpg")

# Section 4: Initialize Parameters ------------------------------------------------------------------------
# * is used for element-wise multiplication while %*% is used for matrix multiplication
radj <- adjacencybuilder(points$id, redgelist)
swadj <- adjacencybuilder(points$id, swedgelist)
beta <- 0.5
epsilon <- 1e-4
theta<- 0.5     # baseline hazard rate
lambda0 <- 1
norm_epsilon <- epsilon / (epsilon + beta)
noedge<- matrix(0, ncol = num_nodes, nrow = num_nodes)
noedge[,as.logical(stats::rbinom(num_nodes, 1, prob = norm_epsilon))] <- 1 # set few nodes as out-sources with the probability norm_epsilon

# Section 5: Integreated Parameter Adjustment -------------------------------------------------------------
lambda1.vec <- c(-5e-5, -3e-5, -2e-5, -1e-5, -9e-6, -8e-6, -7e-6, -6e-6,-5e-6,-4e-6,-3e-6,-2e-6)
maxtime.vec <- c(500, 1000, 2000, 4000, 8000, Inf)
ratemat.lst <- lapply(lambda1.vec, function(x) theta*exp(lambda0+x*dist.mat))
difftimes.lst <- list()                 # set the diagnal of ratemat.lst as 0 before calculate diffusion time
difftimes.df <- rdifftimes.df <- swdifftimes.df <- rg.df <- swg.df <- rsimout.df <- swsimout.df <- 
  data.frame(matrix(vector('list'), length(lambda1.vec), length(maxtime.vec)))
rcas.length <- swcas.length <- data.frame(matrix(vector(), length(lambda1.vec), length(maxtime.vec)))
for (i in 1:length(ratemat.lst)){
  diag(ratemat.lst[[i]]) <- 0
  difftimes.lst[[i]] <- matrix(stats::rexp(num_nodes^2, rate = ratemat.lst[[i]]), nrow = num_nodes)
  rownames(difftimes.lst[[i]]) <- colnames(difftimes.lst[[i]]) <- points$id
  diag(difftimes.lst[[i]]) <- 0
  for (j in 1: length(maxtime.vec)){    # set different cutoff values for the ith lambda1
    tmp <- difftimes.lst[[i]] 
    tmp[tmp >= maxtime.vec[j]] <- 0
    difftimes.df[[i,j]] <- tmp
    rdifftimes.df[[i,j]]<- (radj - noedge)^2 * difftimes.df[[i,j]]
    rg.df[[i,j]] <- igraph::graph.adjacency(rdifftimes.df[[i,j]], weighted = TRUE, mode = "directed")
    swdifftimes.df[[i,j]]<- (swadj - noedge)^2 * difftimes.df[[i,j]]
    swg.df[[i,j]] <- igraph::graph.adjacency(swdifftimes.df[[i,j]], weighted = TRUE, mode = "directed")
    for(k in points$id){                # set each node as initial infected node
      rsimout.df[[i,j]] <- rbind(rsimout.df[[i,j]], simulatecascades(k, rg.df[[i,j]]))
      swsimout.df[[i,j]] <- rbind(swsimout.df[[i,j]], simulatecascades(k, swg.df[[i,j]]))
    }
    rcas.length[i,j] <- nrow(rsimout.df[[i,j]])
    swcas.length[i,j] <- nrow(swsimout.df[[i,j]])
  }
}
rownames(rcas.length) <- rownames(swcas.length) <- lambda1.vec
colnames(rcas.length) <- colnames(swcas.length) <- maxtime.vec
rcas.length.long <- rcas.length %>% rownames_to_column(.,var="lambda1") %>% melt(., id="lambda1", var="maxtime")
swcas.length.long <- swcas.length %>% rownames_to_column(.,var="lambda1") %>% melt(., id="lambda1", var="maxtime")
caslength.df <- data.frame( "lambda1"= as.numeric(rcas.length.long[,1]), "maxtime"=as.numeric(as.character(rcas.length.long[,2])), "rcas.length"= rcas.length.long[,3], "swcas.length" = swcas.length.long[,3])


# Section 6: Visualize Cascades ------------------------------------------------------
save(swsimout.df, file="swsimout.df")


ui = fluidPage(
  titlePanel("Cascades Visualization"),
  fluidRow(
    column(12,
           h4("Geo-visualization of Dynamic Cascades"),
           p(textOutput('dynamicText')),
           leafletOutput("map", height = 650),
           br()
           ),
    column(4,
           numericInput("nodeid", "Initial Infected Node", 1, min=1, max = num_nodes),
           radioButtons("networktype", "Network Structure", c("Random Network","Small-world Network"), selected = NULL, inline=TRUE),
           actionButton("start",label="Start",icon=icon("play","fa-3x"),style="border-color: transparent"),
           div(style="display: inline-block; width: 30px;", HTML("<br>")),
           actionButton("clear",label="Clear", icon=icon("stop","fa-3x"), style = "border-color:transparent")
    ),
    #offset parameter is used on the center column to provide custom spacing
    column(4, selectInput("lambda1", "Choose lambda1",choices = lambda1.vec, selected = NULL, multiple = F)),
    column(4, selectInput("maxtime", "Choose cutoff maximum time",choices = maxtime.vec, selected = NULL, multiple = F))
    )
)
server = function(input, output){
  # A reactive expression uses widget input and returns a value and update the value whenever the original widget changes
  nodesgroup <- as.character(1:num_nodes)
  edgesgroup <- as.character((num_nodes+1):num_edges)
  edges <-  reactive({
    if(input$networktype == "Random Network"){
      rspLines.df
    }else{
      swspLines.df
    }
  })
  cascade <-  reactive({
    lambda1.id <- which(input$lambda1 == lambda1.vec)
    maxtime.id <- which(input$maxtime == maxtime.vec)
    node.id <- input$nodeid
    if(input$networktype == "Random Network"){
      cascades <- rsimout.df[[lambda1.id, maxtime.id]]
      cascade <- cascades[which(cascades$cascade_id == node.id),]
    }else{
      cascades <- swsimout.df[[lambda1.id, maxtime.id]]
      cascade <- cascades[which(cascades$cascade_id == node.id),]
    }
  })
  output$dynamicText <- renderText({
    sprintf('Simulated cascade when the initial infected node is %s, lambda1 is %s, and the cutoff time is %s', input$nodeid, input$lambda1, input$maxtime)
  })
  output$map <- renderLeaflet({
    leaflet() %>% 
      addTiles('http://{s}.basemaps.cartocdn.com/dark_all/{z}/{x}/{y}.png', attribution='Map tiles by <a href="http://stamen.com">Stamen Design</a>, 
               <a href="http://creativecommons.org/licenses/by/3.0">CC BY 3.0</a> &mdash; Map data &copy; <a href="http://www.openstreetmap.org/copyright">OpenStreetMap</a>') %>%
      addCircles(data = sppoints.df, popup=~paste("Node: ", id), radius=10000, stroke=TRUE, color = "red", opacity = 1, weight = 1, fillOpacity = 1, group= nodesgroup) %>%  
      setView(lng = -97.5, lat = 37, zoom = 5)
  })
  observe({
    # Don't miss the open/close brackets after the reactive variable name
    leafletProxy("map") %>% clearGroup(edgesgroup) %>% addPolylines(data= edges(), color = "white", weight = 0.5, group=edgesgroup)
  })
  observeEvent(input$start,{
    for (i in 1:nrow(cascade())){
      infnode <- as.character(cascade()$node_name[i])
      infnode.coor<- sppoints.df@coords[as.numeric(infnode),]
      if(i==1){
        leafletProxy("map") %>% clearGroup(infnode) %>% addCircles(lng=infnode.coor[1], lat=infnode.coor[2], popup=paste("Node: ", infnode), radius=10000, stroke=TRUE, color = "green", opacity = 1, weight = 1, fillOpacity = 1,group=infnode)
      }else{
        parentnode<- cascade()$parent_node[i]
        edgeid <- which(edges()@data[,1]==parentnode & edges()@data[,2]==infnode)
        edgegroup<- as.character(edgeid + num_nodes)
        leafletProxy("map") %>% clearGroup(infnode) %>% clearGroup(edgegroup) %>% 
          addCircles(lng=infnode.coor[1], lat=infnode.coor[2], popup=paste("Node: ", infnode), radius=10000, stroke=TRUE, color = "green", opacity = 1, weight = 1, fillOpacity = 1,group=infnode) %>%
          addPolylines(data= edges()@lines[[edgeid]], color = "yellow", weight = 2,group=edgegroup)
      }
      #if(i !=nrow(cascade())) Sys.sleep((cascade()$event_time[i+1]- cascade()$event_time[i])/1000) # set sleep time
    }
  })
  
  observeEvent(input$clear,{
    # Don't miss the open/close brackets after the reactive variable name
    leafletProxy("map") %>% clearShapes() %>% 
      addCircles(data = sppoints.df, popup=~paste("Node: ", id), radius=10000, stroke=TRUE, color = "red", opacity = 1, weight = 1, fillOpacity = 1, group= nodesgroup) %>%
      addPolylines(data= edges(), color = "white", weight = 0.5,group=edgesgroup)
  })
}

shinyApp(ui, server)

# Section 7: Statistic Plots -------------------------------------------------------------------------------
### Distance histogram of Random and Small-world Network
rdistance <- as.vector(radj*dist.mat)
rdistance <- data.frame("edgedistance"= sort(rdistance[-which(rdistance ==0)]),"networktype" = "random")
swdistance<- as.vector(swadj*dist.mat)
swdistance <- data.frame("edgedistance"= sort(swdistance[-which(swdistance ==0)]),"networktype" = "small-world")
min(rdistance$edgedistance)
edge.dist <- rbind(rdistance, swdistance)

# calculate the average distance of each group
mean.dist <- ddply(edge.dist, "networktype", summarise, mean=mean(edgedistance), sd = sd(edgedistance))
histogram <- ggplot(edge.dist,aes(x=edgedistance/1000,color= networktype,fill= networktype))+ 
  facet_grid(.~networktype)+
  geom_histogram(position="identity", alpha=0.7,binwidth=70) + 
  geom_vline(data = mean.dist, aes(xintercept=mean/1000), color="black", linetype="dashed", size=0.8,alpha=0.7)+
  geom_density(aes(y =  ..density..*(70*600)), color="black",fill="transparent", size=0.8)+
  scale_color_manual("network",values=c("red", "blue"))+
  scale_fill_manual("network",values=c("red", "blue"))+
  labs(x="Edge Distance (km)", y="Frequency") + 
  theme(axis.title.x = element_text(size=10, face="bold"),
        axis.title.y = element_text(size=10, face="bold"), legend.position = c(0.07, 0.9),
        legend.title = element_text(size = 10, face = "bold"), legend.text = element_text(size = 8),
        legend.background = element_rect(fill="transparent"))

qqnorm<- ggplot(edge.dist,aes(sample=edgedistance/1000,color= networktype)) + facet_grid(.~networktype) +
  stat_qq() + stat_qq_line(color="black",size=0.8) +
  scale_color_manual("network",values=c("red", "blue"))+
  labs(x="Quantiles", y="Edge Distance (km)")+
  theme(axis.title.x = element_text(size=10, face="bold"),
        axis.title.y = element_text(size=10, face="bold"), legend.position = c(0.07, 0.9),
        legend.title = element_text(size = 10, face = "bold"), legend.text = element_text(size = 8),
        legend.background = element_rect(fill="transparent"))

jpeg("Distance.jpg", width = 5700, height = 6000, units = "px", quality = 100, res=600)
grid.arrange(histogram, qqnorm, ncol = 1, nrow = 2, top = textGrob("Statistical Result of the Distance", gp=gpar(fontface="bold")))
dev.off()

bin.no <- length(which(ggplot_build(histogram)$data[[1]]$group==1))
# count = density * number of datapoints (600)*binwidth(70)

### Scatter plot for different Lambda1 
lambdaplot<- ggplot(caslength.df)+ geom_point(aes(x=lambda1, y=rcas.length, color="rcolor", shape="rshape", size = "rsize")) + 
  geom_line(aes(x=lambda1, y=rcas.length, color="rcolor")) +
  geom_point(aes(x=lambda1, y=swcas.length, color="swcolor",shape="swshape", size = "swsize")) + 
  geom_line(aes(x=lambda1, y=swcas.length, color="swcolor")) +
  #ggtitle(expression("Simulated Cascades for different "~lambda[1]))+
  labs(x=expression(lambda[1]), y="Total Number of Infected Nodes") + 
  scale_color_manual(name = "Network Structure", values = c("rcolor" = "red", "swcolor"="blue"), labels = c("Random Network", "Small-world Network"))+
  scale_shape_manual(name = "Network Structure", values = c("rshape" = 17, "swshape"= 19), labels = c("Random Network", "Small-world Network")) + 
  scale_size_manual(name = "Network Structure", values = c("rsize" = 1, "swsize"= 1), labels = c("Random Network", "Small-world Network")) +
  theme_gray(base_size = 10) + theme(axis.title.x = element_text(size=12, face="bold"),axis.title.y = element_text(size=12, face="bold"), 
                                     legend.position = c(0.25, 0.94), legend.title = element_text(size = 8, face = "bold"), 
                                     legend.text = element_text(size = 8)) + 
  facet_grid(maxtime ~ .)
  # annotation_custom(tableGrob(caslength.df, rows=NULL, theme = ttheme_default(base_size = 7)),xmin = -4.5e-5, xmax = -3.5e-5, ymin = 5000, ymax = 8700)
maxplot <- ggplot(caslength.df[which(caslength.df$maxtime!= Inf),])+ geom_point(aes(x=maxtime, y=rcas.length, color="rcolor", shape="rshape", size = "rsize")) + 
  geom_line(aes(x=maxtime, y=rcas.length, color="rcolor")) +
  geom_point(aes(x=maxtime, y=swcas.length, color="swcolor",shape="swshape", size = "swsize")) + 
  geom_line(aes(x=maxtime, y=swcas.length, color="swcolor")) +
  #ggtitle("Simulated Cascades for different maximum cutoff values")+
  labs(x="maxtime", y="Total Number of Infected Nodes") + 
  scale_color_manual(name = "Network Structure", values = c("rcolor" = "red", "swcolor"="blue"), labels = c("Random Network", "Small-world Network"))+
  scale_shape_manual(name = "Network Structure", values = c("rshape" = 17, "swshape"= 19), labels = c("Random Network", "Small-world Network")) + 
  scale_size_manual(name = "Network Structure", values = c("rsize" = 2, "swsize"= 2), labels = c("Random Network", "Small-world Network")) +
  theme_gray(base_size = 10) + theme(axis.title.x = element_text(size=12, face="bold"), axis.title.y = element_text(size=12, face="bold"), 
                                     axis.text.x = element_text(angle = 60), legend.position = c(0.1, 0.95), 
                                     legend.title = element_text(size = 8, face = 'bold'), legend.text = element_text(size = 8)) + 
  facet_grid(. ~ lambda1)
# save as .jpg
jpeg("SimulatedCascadeLength.jpg", width = 8000, height = 6000, units = "px", quality = 100, res=600)
grid.arrange(lambdaplot, maxplot, ncol = 2, widths = c(3, 7),top = textGrob("Simulated Cascades for different lambda1 and cutoff values of the diffusion time", gp=gpar(fontface="bold")))
dev.off()
### Box-plot function
theme.boxplot <- theme(axis.text.x = element_text(colour = "grey20", size = 12),
                   axis.text.y = element_text(colour = "grey20", size = 12),
                   axis.title.x = element_text(size=12, face="bold"),
                   axis.title.y = element_text(size=12, face="bold"))
rboxplot.lambda<- ggplot(caslength.df[which(caslength.df$maxtime!= Inf & caslength.df$lambda1!= -7e-5),], mapping = aes(x = lambda1, y = rcas.length, group = lambda1))+
  geom_boxplot(fill = "white", colour = "red", outlier.color= "black",outlier.fill = "black", outlier.size = 1.5, outlier.shape = 16) + 
  stat_boxplot(geom ="errorbar", colour = "red") + geom_point(size = 1, color="black", stroke=0) + 
  labs(x=expression(lambda[1]), y="Total Number of Infected Nodes") + 
  theme.boxplot + annotate("text", x = -4e-05, y = 9500, label = "Random Network", size =5, colour = "red", fontface = "bold")
rboxplot.max<- ggplot(caslength.df[which(caslength.df$maxtime!= Inf & caslength.df$lambda1!= -7e-5),], mapping = aes(x = maxtime, y = rcas.length, group = maxtime))+
  geom_boxplot(fill = "white", colour = "red", outlier.color= "black",outlier.fill = "black", outlier.size = 1.5, outlier.shape = 16) + 
  stat_boxplot(geom ="errorbar", colour = "red") + geom_point(size = 1, color="black", stroke=0) + 
  labs(x=expression(maxtime), y=NULL) + theme.boxplot 
swboxplot.lambda<- ggplot(caslength.df[which(caslength.df$maxtime!= Inf & caslength.df$lambda1!= -7e-5),], mapping = aes(x = lambda1, y = swcas.length, group = lambda1))+
  geom_boxplot(fill = "white", colour = "blue", outlier.color= "black",outlier.fill = "black", outlier.size = 1.5, outlier.shape = 16) + 
  stat_boxplot(geom ="errorbar", colour = "blue") + geom_point(size = 1, color="black", stroke=0) +
  labs(x=expression(lambda[1]), y="Total Number of Infected Nodes") + 
  theme.boxplot + annotate("text", x = -3.8e-05, y = 9000, label = "Small-world Network", size =5, colour = "blue", fontface = "bold")
swboxplot.max<- ggplot(caslength.df[which(caslength.df$maxtime!= Inf & caslength.df$lambda1!= -7e-5),], mapping = aes(x = maxtime, y = swcas.length, group = maxtime))+
  geom_boxplot(fill = "white", colour = "blue", outlier.color= "black",outlier.fill = "black", outlier.size = 1.5, outlier.shape = 16) + 
  stat_boxplot(geom ="errorbar", colour = "blue") + geom_point(size = 1, color="black", stroke=0) +
  labs(x=expression(maxtime),y= NULL) + theme.boxplot 
jpeg("Boxplots.jpg", width = 6000, height = 6000, units = "px", quality = 100, res=600)
grid.arrange(rboxplot.lambda, rboxplot.max, swboxplot.lambda, swboxplot.max, ncol = 2, nrow=2,
             widths = c(5, 5), top = textGrob("Boxplots of the number of infected nodes for different lambda1 and cutoff values of the diffusion time",gp=gpar(fontface="bold")))
dev.off()

# # Infer the network -------------------------------------------------------
# maxtime <- 2000
# rlambda1<- -2e-05
# swlambda1<- -9e-06
# try <- rsimout.df[[5,3]]
# swcascade <-  swsimout.df[[3,3]]
# rcascade.l <- as_cascade_long(rcascade)
# swcascade.l <- as_cascade_long(swcascade)
# 
# # Generate by Rcpp::compileAttributes()
# netinf_ <- function(cascade_nodes, cascade_times, n_edges, model, params, quiet, auto_edges, cutoff) {
#   .Call(`_NetworkInference_netinf_`, cascade_nodes, cascade_times, n_edges, model, params, quiet, auto_edges, cutoff)
# }
# 
# count_possible_edges_ <- function(cascade_nodes, cascade_times, quiet = TRUE) {
#   .Call(`_NetworkInference_count_possible_edges_`, cascade_nodes, cascade_times, quiet)
# }
# 
# # Simulate Cascades -------------------------------------------------------
# coord <- as.matrix(points.df[,-1])
# DIST <- exponentialHaz()
# COVMODEL <- ExponentialCovFct() # spatial corvaince function Y_i
# COVPARS <- c(0.7, 0.1) # \sigma and \phi in spatial covariance function
# # for exponential hazard: \omega = \theta in my thesis; for weibull: \omega =(\alpha, \lambda)
# # the \beta here is associated with \lambda in my thesis
# # cbind or rbind can do for only one vector
# sim <- simsurv(X = cbind(constant = rep(1, num_nodes)),
#                beta = 0.0296, dist = DIST, omega = 1,
#                cov.parameters = COVPARS, cov.model = COVMODEL, coord = coord,
#                mcmc.control = mcmcpars(nits = 110000, burn = 10000, thin = 100))
# sim.df <- cbind(sim$survtimes, sim$coords)
# colnames(sim.df)<- c("survtimes","x","y")
# # True simulated survival Time
# survtimes <- sim$survtimes
# # The exponential distribution with rate ?? has density f(x) = ?? {e}^{- ?? x}
# # rexp: random generation for the exponential distribution with rate ??
# # rate = 1/mean(survtimes[,1])
# censtimes <- rexp(num_nodes, 1/mean(survtimes))
# # Right-censored data with survival time and status. If the censoring time is less than 
# # the simulated survival time, the status is 0(alive), otherwise 1(dead)
# simsurv <- gencens(survtimes, censtimes)
# simsurv.df <- cbind(data.frame(as.matrix(simsurv)), sim$X)
# spatcas <- SpatialPointsDataFrame(coord, simsurv.df)
# m <- leaflet() %>% addTiles('http://{s}.basemaps.cartocdn.com/dark_all/{z}/{x}/{y}.png', 
#                             attribution='Map tiles by <a href="http://stamen.com">Stamen Design</a>, 
#                             <a href="http://creativecommons.org/licenses/by/3.0">CC BY 3.0</a> &mdash; Map data &copy; <a href="http://www.openstreetmap.org/copyright">OpenStreetMap</a>') %>%
#   addCircles(data = spatcas, radius=10000, stroke=TRUE, color = "red", opacity = 1, weight = 1, fillOpacity = 1)
# 
# estimate <- survspat(Surv(time, status) ~ age + sex, dist=DIST,
#                      cov.model = COVMODEL,data=spatdat,
#                      mcmc.control = mcmcpars(nits = 500000, burn = 10000, thin = 490))
# # Define Edge Score
# edgeScore <- function (cascade_times, pairid){
#   timeinterval <- cascade_times[pairid[2]]- cascade_times[pairid[1]]
# }
# # Find possible edges for all cascades
# getpossibleEdges <- function (cascade){
#   n_cascade <- length(cascade$cascade_nodes)
#   possibleedges.lst <- list(data.frame(matrix(ncol=2,nrow = 0)))
#   tmp <- data.frame(matrix(ncol=2,nrow = 0))
#   for(i in 1:n_cascade){
#     nsize <- length(cascade$cascade_nodes[[i]]) # node size in ith cascade
#     for (j in 1:(nsize-1)){
#       parent <- cascade$cascade_nodes[[i]][j]
#       t_parent <- cascade$cascade_times[[i]][j]
#       for(k in (j+1): nsize){
#         child <- cascade$cascade_nodes[[i]][k]
#         t_child <- cascade$cascade_times[[i]][k]
#         if(t_parent == t_child) next
#         pair_id <- c(parent, child)
#         score <- edgeScore(pair_id)
#         tmp <- rbind(tmp, pair_id, stringsAsFactors = FALSE)
#       }
#     }
#     colnames(try) <- c("parent","child")
#     possibleedges.lst[[i]] <- try
#   }
#   return(possibleedges.lst)
# }
### reset the colname values in the nested df
# for (i in 1:12){
#   for(j in 1:6){
#     colnames(swsimout.df[[i,j]])[colnames(swsimout.df[[i,j]])=="parentnode"] <- "parent_node"
#   }
# }
# 
# # Diagnose codes for parameters
# diag(ratemat.lst[[1]]) <- 0
# difftimes.lst[[1]] <- matrix(stats::rexp(num_nodes^2, rate = ratemat.lst[[1]]), nrow = num_nodes)
# rownames(difftimes.lst[[1]]) <- colnames(difftimes.lst[[1]]) <- points$id
# diag(difftimes.lst[[1]]) <- 0
# tmp <- difftimes.lst[[1]] 
# tmp[tmp >= maxtime.vec[1]] <- 0
# difftimes.df[[1,1]] <- tmp
# rdifftimes.df[[1,1]]<- (radj - noedge)^2 * difftimes.df[[1,1]]
# rg.df[[1,1]] <- igraph::graph.adjacency(rdifftimes.df[[1,1]], weighted = TRUE, mode = "directed")
# shortestpath <- as.vector(igraph::shortest_paths(rg.df[[1,1]], from = 1, mode = "out", output = "vpath", predecessors=TRUE))
# parentnode <- as.vector(shortestpath$predecessors)
# dists <- as.vector(igraph::distances(rg.df[[1,1]], v = 1, mode = "out"))
# tmp <- data.frame("node_name" = V(rg.df[[1,1]])$name, "event_time"= as.numeric(dists),"parent_node" = parentnode, "cascade_id" = 1,stringsAsFactors=FALSE)
# tmp <- tmp[-which(tmp$event_time==Inf),]
# rsimout.df[[1,1]] <- rbind(rsimout.df[[1,1]], simulatecascades(k, rg.df[[1,1]]))
