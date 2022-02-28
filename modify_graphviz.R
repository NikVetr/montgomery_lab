# set our own layout using ggraph
# see http://blog.schochastics.net/post/ggraph-tricks-for-common-problems/

# when using source, this will fail if one of the libraries was not installed
library(ggraph)
library(scatterpie)
library(grid)
library(networkD3)
library(igraph)
library(repfdr)
library(data.table)

#' An auxiliary function for reducing a list of sets by a regex.
#' Useful for filtering DEA analyte sets by tissues or omes.
#' 
#' @params l a named list of character vectors
#' @params regs a character vector of regular expressions
limit_sets_by_regex<-function(sets,regs){
  if(is.null(regs) || length(regs)==0){return(sets)}
  l = list()
  for(nn in names(sets)){
    v = sets[[nn]]
    newv = c()
    for(r in regs){
      r = paste0(r,";")
      newv = union(newv,v[grepl(r,v)])
    }
    l[[nn]] = newv
  }
  return(l)
}

#' Auxiliary function to get the largest paths in a graph.
#' 
#' This is implemented using a dynamic programming approach where
#' we iteratively add the data of the next edge.
#' 
#' When analyzing an edge (x,y) with a set of analytes s, we extend all
#' trajectories that end with x using the new edge, but also all the 
#' trajectories that start with y. 
#' 
#' At the end, because we examine all edges and all extensions we are
#' guaranteed to have covered all full paths.
#' 
#' The min_size parameter is important, since we are interested in paths
#' and edges with at least this number of analytes, then we know that if a
#' current trajectory does not have at least min_size analytes, then since the 
#' set of any extension can only be the same or smaller then we can ignore
#' such paths moving forward. 
#' 
#' However, this parameter has to be considered with care, as specifying a number
#' that is too high will result in no paths in the output.
#' 
#' @param edge_sets a named list of string vectors. The name of an edge is node_id---node_id
#'        edges with no analytes have a NULL set (a set of size zero, but are still represented),
#'        node ids are time_points_Fx_My where x and y represent the up/down state in each sex.
#' @param min_size a number, specifying the minimal path size to be considered
#' @return NULL if no paths of size of at least min_size were found, otherwise
#'         return a data frame that represents all paths of size min_size or greater,
#'         ranked from the largest path to the smallest one.
get_trajectory_sizes_from_edge_sets<-function(edge_sets,min_size=10){
  node_names = unique(unlist(strsplit(names(edge_sets),split="---")))
  node_weeks = sapply(node_names,function(x)strsplit(x,split="_")[[1]][1])
  full_path_size = length(unique(node_weeks))
  arrs = strsplit(names(edge_sets),split="---")
  prevs = sapply(arrs,function(x)x[1])
  nexts = sapply(arrs,function(x)x[2])
  l = list()
  l_sets = c()
  for(j in 1:length(edge_sets)){
    if(length(edge_sets[[j]])<min_size){next}
    ind = length(l)+1
    curr_prev = prevs[j]
    curr_next = nexts[j]
    l[[ind]] = c(curr_prev,curr_next)
    l_sets[[ind]] = edge_sets[[j]]
    for(j2 in 1:length(l)){
      j2_arr = l[[j2]]
      # here the path of j2 starts with the current edge next
      # so the edge represented by j can extend the path j2
      # "on the left" as it represent a previous time interval
      if(j2_arr[1] == curr_next){
        new_ind = length(l)+1
        l[[new_ind]] = c(curr_prev,j2_arr)
        l_sets[[new_ind]] = intersect(l_sets[[j2]],edge_sets[[j]])
      }
      # here the path of j2 end with the current edge prev/start
      # so the edge represented by j can extend the path j2
      # "on the right" as it represent a subsequent time interval
      if(j2_arr[length(j2_arr)]==curr_prev){
        new_ind = length(l)+1
        l[[new_ind]] = c(j2_arr,curr_next)
        l_sets[[new_ind]] = intersect(l_sets[[j2]],edge_sets[[j]])
      }
    }
    keep = sapply(l_sets,length) >= min_size
    l = l[keep]
    l_sets = l_sets[keep]
  }
  # limit the results to the "full" trajectories
  traj_sizes = sapply(l,length)
  # here we assume that at least one full pathway survived
  # the min size filter above, otherwise, return NULL
  traj_inds = traj_sizes==full_path_size
  if(sum(traj_inds)==0){return(NULL)}
  l = l[traj_inds]
  l_sets = l_sets[traj_inds]
  trajectories = cbind(
    t(sapply(l,function(x)x)),
    sapply(l_sets,length)
  )
  trajectories = data.frame(trajectories,stringsAsFactors = F)
  trajectories[[6]] = as.numeric(trajectories[[6]])
  trajectories = trajectories[order(trajectories[[6]],decreasing = T),]
  return(trajectories)
}

#' Keep the edges of the top trajectories of an edge set
filter_edge_sets_by_trajectories<-function(edge_sets,topk=5,min_path_size=5){
  
  traj = get_trajectory_sizes_from_edge_sets(edge_sets,min_size = min_path_size)
  edges_to_keep = c()
  if(!is.null(traj)){
    topk = min(topk,nrow(traj))
    traj = traj[1:topk,]
    for(j in 1:(ncol(traj)-1)){
      curr_edges = paste(traj[,j],traj[,j+1],sep="---")
      edges_to_keep = union(edges_to_keep,curr_edges)
    }
  }
  
  e_copy = copy(edge_sets)
  for(e in names(e_copy)){
    if(!(e %in% edges_to_keep)){
      e_copy[[e]] = character(0) 
    }
  }
  return(e_copy)
}

#' The main function for obtaining a graphical (tree) representation of the differential
#' analysis results.
#' 
#' @param tissues acharacter vector. The set of tissues to take for the analysis. If null take all.
#' @param omes acharacter vector. The set of omes to take for the analysis. If null take all.
#' @param node_sets a named list with the node (state) sets of analyte, see details for analyte name convention.
#' @param edge_sets a named list with the edge (state) sets of analyte, see details for analyte name convention.
#' @param min_size a numeric, a threshold on the set sizes to be considered
#' @param parallel_edges_by_ome a logical. TRUE means that we want to added parallel edges for the different omes.
#' @param parallel_edges_by_tissue a logical. TRUE means that we want to added parallel edges for the different tissues.
#' @param edge_width_range a numeric vector of size 2, a parameter for ggraph
#' @param edge_alpha_range a numeric vector of size 2, a parameter for ggraph
#' @param color_nodes_by_states a logical. If TRUE, nodes are colored by states. 
#'        Red for up-reg, blue for down-reg, green for a discrepancy between the sexes.
#' @param max_trajectories a numeric or NULL. If not NULL then it specifies the number of pathways
#'        to keep when looking into the edge sets after filtering by omes and tissues. If 
#'        parallel_edges_by_tissue = TRUE then take the top trajectories in each tissue.
#' @param highlight_subset a character string or NULL. If not NULL then it specifies the name of a 
#'        node, edge, or path to highlight in the tree. 
#' 
#' 
#' @details 
#'         The function filters the input set to include analytes from the given tissues
#'         and omes (if tissues/omes are not null).
#'         If parallel edges are requested then the relevant edge sizes are computed internally
#'         and are used within ggraph for the output plot.
#'         Analyte names are in the ome;tissue;feature_id format.
#'         Node set names are are in the week_number(1,2,4,8)'w'_F{-1,0,1}_M{-1,0,1} format
#'         for example: 1w_F-1_M0
#'         Edge set names are in the node_a---node_b format, e.g., 4w_F0_M0---8w_F0_M1 
get_tree_plot_for_tissue<-function(
  tissues,
  omes = NULL,
  node_sets,
  edge_sets,
  min_size = 20,
  parallel_edges_by_ome = F,
  parallel_edges_by_tissue = F,
  edge_width_range = c(0,10),
  edge_alpha_range = c(0,1),
  color_nodes_by_states = T,
  max_trajectories = NULL,
  highlight_subset = NULL
){
  
  tissues = unique(tissues)
  omes = unique(omes)
  
  tissue_node_sets = limit_sets_by_regex(node_sets,tissues)
  tissue_node_sets = limit_sets_by_regex(tissue_node_sets,omes)
  tissue_edge_sets = limit_sets_by_regex(edge_sets,tissues)
  tissue_edge_sets = limit_sets_by_regex(tissue_edge_sets,omes)
  
  # Set things up to highlight a subset of the tree
  highlight_edge = F
  highlight_node = F
  highlight_path = F
  if(!is.null(highlight_subset)){
    if(grepl(":", highlight_subset)){
      # Check if it's this tissue 
      if(grepl(":", highlight_subset)){
        if(length(tissues)>1){
          warning("'highlight_subset' cannot handle a subset of tissues. Setting 'highlight_subset' to NULL Specify a single tissue or remove the tissue prefix from 'highlight_subset'.")
          highlight_subset = NULL
        }else{
          highlight_tissue = gsub(":.*","",highlight_subset)
          highlight_cluster = gsub(".*:","",highlight_subset)
          if(highlight_tissue != tissues){
            warning("'highlight_subset' tissue prefix does not match the tissue supplied. Setting 'highlight_subset' to NULL.")
            highlight_subset = NULL
          }
        }
      }else{
        highlight_cluster = highlight_subset
      }
    }
    if(!is.null(highlight_subset)){
      if(grepl("---", highlight_subset)){
        highlight_edge = T
      }else if(grepl("->",highlight_subset)){
        highlight_path = T
      }else{
        highlight_node = T
      }
    }
  }
  
  if(!is.null(highlight_subset) & any(parallel_edges_by_ome, parallel_edges_by_tissue, color_nodes_by_states)){
    warning("Setting 'parallel_edges_by_ome', 'parallel_edges_by_tissue', and 'color_nodes_by_states' to FALSE because 'highlight_subset' is TRUE.")
    parallel_edges_by_ome = F
    parallel_edges_by_tissue = F
    color_nodes_by_states = F
  }
  
  # Check the input omes and tissues sets
  # If either ome or tissues is NULL, take them from the data
  all_analytes = unique(unlist(node_sets))
  all_analytes_info = strsplit(all_analytes,split=";")
  data_tissue_set = unique(sapply(all_analytes_info,function(x)x[2]))
  data_ome_set = unique(sapply(all_analytes_info,function(x)x[1]))
  if(!all(tissues %in% data_tissue_set)){
    warning("Some input tissue names are not in the dataset and will be removed")
    tissues = intersect(tissues,data_tissue_set)
  }
  if(!all(omes %in% data_ome_set)){
    warning("Some input ome names are not in the dataset and will be removed")
    omes = intersect(omes,data_ome_set)
  }
  if(is.null(tissues)){
    tissues = data_tissue_set
  }
  if(is.null(omes)){
    omes = data_ome_set
  }
  
  # if we are not asked to get parallel edges and are asked to take 
  # the top trajectories then we can simply filter the edge sets and move on
  if(!is.null(max_trajectories) && !parallel_edges_by_tissue){
    tissue_edge_sets = filter_edge_sets_by_trajectories(tissue_edge_sets,max_trajectories)
  }
  
  # Create the initial matrices/data frame
  # with some additional info (depends on the user input), these will be used as input
  # for generating the graph
  node_info = cbind(
    names(tissue_node_sets),
    sapply(names(tissue_node_sets),function(x)strsplit(x,split="_")[[1]][1]),
    sapply(tissue_node_sets,length)
  )
  colnames(node_info) = c("node","week","size")
  edge_info = cbind(
    sapply(names(tissue_edge_sets),function(x)strsplit(x,split="---")[[1]][1]),
    sapply(names(tissue_edge_sets),function(x)strsplit(x,split="---")[[1]][2]),
    sapply(tissue_edge_sets,length)
  )
  colnames(edge_info) = c("prev","next","size")
  
  # If we are asked to add parallel edges by omes then we need to 
  # compute ome edge sizes and add them
  if(parallel_edges_by_ome){
    for(o in omes){
      v = c()
      for(e in rownames(edge_info)){
        v[e] = sum(grepl(paste0(o,";"),tissue_edge_sets[[e]]))
      }
      edge_info = cbind(edge_info,v)
      colnames(edge_info)[ncol(edge_info)] = o
    }
    # always use the same assay colors
    edge_colors = define_assay_cols() # defined in cluster_viz_fx.R
    # require(RColorBrewer)
    # o = c("METAB","TRNSCRPT","PROT","ACETYL","PHOSPHO","UBIQ","ATAC","METHYL","IMMUNO")
    # edge_colors = brewer.pal(length(o), 'Set1')
    # names(edge_colors) = o
  }
  # if parallel_edges_by_tissue = TRUE then:
  # if we do not need to filter the top trajectories then we can simply count
  # the number of features per edge for a tissue using regular expressions
  if(is.null(max_trajectories) && parallel_edges_by_tissue){
    for(tissue in tissues){
      v = c()
      for(e in rownames(edge_info)){
        v[e] = sum(grepl(paste0(tissue,";"),tissue_edge_sets[[e]]))
      }
      edge_info = cbind(edge_info,v)
      colnames(edge_info)[ncol(edge_info)] = tissue
    }
    # always use the same assay colors
    edge_colors = MotrpacBicQC::tissue_cols
  }
  # if parallel_edges_by_tissue = TRUE then:
  # if we are asked to filter by the top trajectories then we need to 
  # recompute the edge sizes per tissue by limiting the tissue edge sets
  if(!is.null(max_trajectories) && parallel_edges_by_tissue){
    for(tissue in tissues){
      curr_tissue_edge_sets = limit_sets_by_regex(tissue_edge_sets,tissue)
      curr_tissue_edge_sets = filter_edge_sets_by_trajectories(
        curr_tissue_edge_sets,max_trajectories)
      v = c()
      for(e in rownames(edge_info)){
        v[e] = length(curr_tissue_edge_sets[[e]])
      }
      edge_info = cbind(edge_info,v)
      colnames(edge_info)[ncol(edge_info)] = tissue
    }
    # always use the same assay colors
    edge_colors = MotrpacBicQC::tissue_cols
  }
  
  # transform the initial node/edge info matrices to data frames
  # filter edge sizes by the specified user input
  d = as.data.frame(edge_info)
  for(j in 3:ncol(d)){
    d[[j]] = as.numeric(d[[j]])
    d[d[,j]<min_size,j] = 0 # filter by min edge size
  }
  d_nodes = data.frame(node_info,check.names = F)
  d_nodes[[3]] = as.numeric(d_nodes[[3]])
  d_nodes$inds = 0:(nrow(d_nodes)-1)
  
  # If we need to plot parallel edges then we need to transform the edge-based
  # data frame into a long format (to make it work with ggraph)
  if(parallel_edges_by_ome || parallel_edges_by_tissue){
    added_cols = names(d)[-c(1:3)]
    if(!parallel_edges_by_ome){added_cols = setdiff(added_cols,omes)}
    if(!parallel_edges_by_tissue){added_cols = setdiff(added_cols,tissues)}
    newd = d[,1:2]
    newd$size = 0
    newd$type = NA
    for(o in added_cols){
      currd = d[,c(1:2,which(colnames(d)==o))]
      names(currd)[3] = "size"
      currd$type = o
      newd = rbind(newd,currd)
    }
    newd = newd[!is.na(newd$type),]
    d_g = igraph::graph_from_data_frame(newd)
    if(all(E(d_g)$size == 0)){
      message("No non-0 edges. Skipping.")
      return()
    }
    E(d_g)$edge_size = edge_attr(d_g,"size")
    d_g_auto_layout <- create_layout(d_g, layout = 'auto')
  }else{
    d_g = igraph::graph_from_data_frame(d)
    if(all(E(d_g)$size == 0)){
      message("No non-0 edges. Skipping.")
      return()
    }
    E(d_g)$edge_size = edge_attr(d_g,"size")
    d_g_auto_layout <- create_layout(d_g, layout = 'auto')
  }
  
  ################
  # From this point we start generating ggraph info
  # (e.g., layout, labels etc.)
  d_g_ordered_nodes = c("F-1_M1","F1_M-1","F0_M-1","F-1_M-1","F-1_M0",
                        "F0_M0","F0_M1","F1_M1","F1_M0")
  d_g_ordered_cols = c("lightgreen","lightgreen","blue","blue","blue",
                       "gray","red3","red3","red3")
  #d_g_ordered_cols_alt = c("black","black","black","black","black",
  #                     "gray","black","black","black")
  d_g_ordered_cols_alt = c("gray4","gray4","gray4","gray4","gray4",
                           "gray","gray4","gray4","gray4")
  d_g_ordered_shapes = c(18,18,19,19,19,19,19,19,19)
  layer_plot_names = c(
    "F down, M up","F up, M down",
    "M only down","Both down","F only down",
    "No\nchange",
    "M only up","Both up","F only up"
  )
  
  names(d_g_ordered_cols) = d_g_ordered_nodes
  names(d_g_ordered_shapes) = d_g_ordered_nodes
  names(d_g_ordered_cols_alt) = d_g_ordered_nodes
  l_x_lim = c(min(d_g_auto_layout$x),max(d_g_auto_layout$x))
  l_y_lim = c(min(d_g_auto_layout$y),max(d_g_auto_layout$y))
  d_g_our_layout  = copy(d_g_auto_layout)
  # set 5 horiz layers over time
  xjump = (l_x_lim[2]-l_x_lim[1]) / 4
  # set 9 vertical layers over ordered_nodes
  yjump = (l_y_lim[2]-l_y_lim[1]) / 8
  d_g_our_layout[d_g_our_layout$name=="0w","x"] = l_x_lim[1]
  d_g_our_layout[d_g_our_layout$name=="0w","y"] = mean(l_y_lim)
  curr_weeks = c("1w","2w","4w","8w")
  for(j in 1:length(curr_weeks)){
    w = curr_weeks[j]
    d_g_our_layout[grepl(w,d_g_our_layout$name),"x"] = l_x_lim[1] + j*xjump
  }
  for(j in 1:length(d_g_ordered_nodes)){
    n = d_g_ordered_nodes[j]
    d_g_our_layout[grepl(n,d_g_our_layout$name),"y"] = l_y_lim[1] + (j-1)*yjump
  }
  # make sure that the 0w node is in the same line as of the no response
  d_g_our_layout[1,"y"] = d_g_our_layout[grepl("F0_M0",d_g_our_layout$name),"y"][1]
  # set node sizes and other features
  V(d_g)$setsize = d_nodes[V(d_g)$name,"size"]
  V(d_g)$setsize[V(d_g)$name == "0w"] = median(V(d_g)$setsize)
  V(d_g)$label = sapply(V(d_g)$name,
                        function(x){a=strsplit(x,split="w_")[[1]];a[length(a)]})
  V(d_g)$label = gsub("_","\n",V(d_g)$label)
  V(d_g)$col = "gray"
  V(d_g)$alt_col = "gray"
  V(d_g)$shape = 15
  for(j in 1:length(d_g_ordered_nodes)){
    n = d_g_ordered_nodes[j]
    V(d_g)$col[grepl(n,d_g_our_layout$name)] = d_g_ordered_cols[n]
    V(d_g)$alt_col[grepl(n,d_g_our_layout$name)] = d_g_ordered_cols_alt[n]
    V(d_g)$shape[grepl(n,d_g_our_layout$name)] = d_g_ordered_shapes[n]
  }
  # for(ome in colnames(tissue_ome_data)){
  #   d_g = set_vertex_attr(d_g,ome,V(d_g),tissue_ome_data[V(d_g)$name,ome])
  # }
  grid_group_annotation_y = c()
  for(j in 1:length(d_g_ordered_nodes)){
    n = d_g_ordered_nodes[j]
    grid_group_annotation_y[n] = l_y_lim[1] + (j-1)*yjump
  }
  weeks_annotation_x = c()
  for(j in 1:length(curr_weeks)){
    w = curr_weeks[j]
    weeks_annotation_x[w] = l_x_lim[1] + j*xjump
  }
  
  # Set colors to highlight edge, node, path
  # Red to highlight; gray for everything else 
  if(highlight_edge){
    V(d_g)$col = "gray"
    V(d_g)$alt_col = "gray"
    # Make edge red
    ecols = rep("gray", nrow(d))
    ecols[which(rownames(d)==highlight_cluster)] = "red"
    E(d_g)$col = ecols
  }else if(highlight_node){
    E(d_g)$col = "gray"
    vcols = rep("gray", length(names(V(d_g))))
    vcols[which(names(V(d_g))==highlight_cluster)] = "red"
    V(d_g)$col = vcols
    V(d_g)$alt_col = vcols
  }else if(highlight_path){
    # extract nodes
    curr_nodes = c("0w",unname(unlist(strsplit(highlight_cluster, "->"))))
    # extract edges
    edge_vector = c()
    for(e in 2:length(curr_nodes)){
      edge_vector = c(edge_vector, sprintf("%s---%s", curr_nodes[e-1], curr_nodes[e]))
    }
    # set colors for edges
    ecols = rep("gray", nrow(d))
    ecols[rownames(d)%in%edge_vector] = "red"
    E(d_g)$col = ecols
    # set colors for nodes
    vcols = rep("gray", length(names(V(d_g))))
    vcols[names(V(d_g))%in%curr_nodes] = "red"
    V(d_g)$col = vcols
    V(d_g)$alt_col = vcols
  }
  
  if(parallel_edges_by_ome || parallel_edges_by_tissue){
    p = ggraph(d_g,layout = d_g_our_layout,) +  
      geom_edge_fan0(
        aes(width = E(d_g)$edge_size,alpha=E(d_g)$edge_size,colour = E(d_g)$type),
        end_cap = circle(3, 'mm'),curvature = 0.1) +
      scale_edge_color_discrete(name="")
  }else if(!is.null(highlight_subset)){
    p = ggraph(d_g,layout = d_g_our_layout,) +  
      #geom_edge_link(aes(width = E(d_g)$edge_size,alpha=E(d_g)$edge_size)) + 
      geom_edge_arc(aes(width = E(d_g)$edge_size,alpha=E(d_g)$edge_size, color = E(d_g)$col),
                    end_cap = circle(3, 'mm'),curvature = 0.1) +
      scale_edge_colour_identity(guide = "none")
  }else{
    p = ggraph(d_g,layout = d_g_our_layout,) +  
      #geom_edge_link(aes(width = E(d_g)$edge_size,alpha=E(d_g)$edge_size)) + 
      geom_edge_arc(aes(width = E(d_g)$edge_size,alpha=E(d_g)$edge_size),
                    end_cap = circle(3, 'mm'),curvature = 0.1) +
      scale_edge_color_discrete(name="")
  }
  
  p = p +
    scale_edge_width(range=edge_width_range,name="Intersect size") + 
    scale_edge_alpha(range=edge_alpha_range,name="Intersect size") 
  
  if(color_nodes_by_states){
    p = p + geom_node_point(aes(size = V(d_g)$setsize),alpha=1,
                            color=V(d_g)$col,shape=V(d_g)$shape)
  }else{
    p = p + geom_node_point(aes(size = V(d_g)$setsize),alpha=1,
                            color=V(d_g)$alt_col,shape=V(d_g)$shape)
  }
  p = p +
    #geom_node_text(aes(label = V(d_g)$label), repel=F) + 
    scale_size(range = c(2,20),name="Number of analytes") +
    theme(panel.background = element_rect(fill="white"),
          legend.key.size = unit(0.4, 'cm'))
  
  names(grid_group_annotation_y) = layer_plot_names
  for(g in names(grid_group_annotation_y)){
    p = p + annotate(geom="text",x=l_x_lim[1] + 0.2, y=grid_group_annotation_y[g],
                     label=g,fontface = "bold")
  }
  names(weeks_annotation_x) = paste("week",c(1,2,4,8))
  for(w in names(weeks_annotation_x)){
    p = p + annotate(geom="text",y=l_y_lim[2]+0.5, x=weeks_annotation_x[w],
                     label=w,fontface = "bold")
  }
  #plot(p)
  
  # change colors and legend edge size
  if(parallel_edges_by_ome | parallel_edges_by_tissue){
    p = p +
      scale_edge_color_manual(values=edge_colors, name="") +
      guides(edge_color=guide_legend(override.aes = list(edge_width=3)))
  }
  
  return(p)
}


#' An auxiliary function for getting the top analyte sets for a set of tissues
extract_tissue_sets<-function(tissues,node_sets,edge_sets,k=3,
                              min_size=20,add_week8=T,omes=NULL){
  
  tissue_node_sets = limit_sets_by_regex(node_sets,tissues)
  tissue_node_sets = limit_sets_by_regex(tissue_node_sets,omes)
  tissue_edge_sets = limit_sets_by_regex(edge_sets,tissues)
  tissue_edge_sets = limit_sets_by_regex(tissue_edge_sets,omes)
  l = list()
  
  # node sets
  node_set_sizes = sapply(tissue_node_sets,length)
  selected_node_sets = names(sort(node_set_sizes,decreasing = T))[1:(k+1)]
  if(add_week8){
    selected_node_sets = union(selected_node_sets,
                               names(tissue_node_sets)[grepl("8w",names(tissue_node_sets))])
  }
  selected_node_sets = setdiff(selected_node_sets,"0w")
  l[selected_node_sets] = tissue_node_sets[selected_node_sets]
  
  # edge sets
  edge_set_sizes = sapply(tissue_edge_sets,length)
  edge_set_sizes = edge_set_sizes[!grepl("0w",names(edge_set_sizes))]
  selected_edge_sets = names(sort(edge_set_sizes,decreasing = T))[1:k]
  l[selected_edge_sets] = tissue_edge_sets[selected_edge_sets]
  
  # path sets: get the top trajectories first
  tissue_top_trajs = get_trajectory_sizes_from_edge_sets(tissue_edge_sets,min_size = min_size)
  if(!is.null(tissue_top_trajs) && nrow(tissue_top_trajs)>0){
    k_paths = min(k,nrow(tissue_top_trajs))
    tissue_top_trajs = tissue_top_trajs[1:k_paths,]
    for(path_i in 1:nrow(tissue_top_trajs)){
      curr_set = tissue_node_sets[[tissue_top_trajs[path_i,2]]]
      for(path_j in 3:5){
        curr_set = intersect(curr_set,
                             tissue_node_sets[[tissue_top_trajs[path_i,path_j]]])
      }
      if(length(curr_set) != tissue_top_trajs[path_i,6]){
        warning("Computed trajectory set size is different from the precomputed size, 
                this happns if edge sets are not simple intersections of node sets")
      }
      l[[paste(tissue_top_trajs[path_i,2:5],collapse="->")]] = curr_set
    }
  }
  l = l[sapply(l,length)>min_size]
  return(l)
}


#' Pull out all non-empty trajectories 
#' 
#' @param node_sets use `load_graph_vis_data()$node_sets`
#' @param edge_sets use `load_graph_vis_data()$edge_sets`
#' @param tissues string vector, optional. tissue subset. all tissues by default
#' @param omes string vector, optional. ome subset. all omes by default
#' 
#' @return named list with one element per trajectories. members are features in the path 
get_all_trajectories = function(edge_sets, 
                                node_sets, 
                                tissues = unique(unname(MotrpacBicQC::tissue_abbr)),
                                omes = unique(unname(MotrpacBicQC::assay_abbr))){
  require("MotrpacBicQC")
  
  l = list()
  tissue_node_sets = limit_sets_by_regex(node_sets,tissues)
  tissue_node_sets = limit_sets_by_regex(tissue_node_sets,omes)
  tissue_edge_sets = limit_sets_by_regex(edge_sets,tissues)
  tissue_edge_sets = limit_sets_by_regex(tissue_edge_sets,omes)
  tissue_top_trajs = get_trajectory_sizes_from_edge_sets(tissue_edge_sets,min_size = 1)
  if(length(tissue_top_trajs) == 0){
    return()
  }
  for(path_i in 1:nrow(tissue_top_trajs)){
    curr_set = tissue_node_sets[[tissue_top_trajs[path_i,2]]]
    for(path_j in 3:5){
      curr_set = intersect(curr_set,
                           tissue_node_sets[[tissue_top_trajs[path_i,path_j]]])
    }
    if(length(curr_set) != tissue_top_trajs[path_i,6]){
      warning("Computed trajectory set size is different from the precomputed size, 
                  this happns if edge sets are not simple intersections of node sets")
    }
    l[[paste(tissue_top_trajs[path_i,2:5],collapse="->")]] = curr_set
  }
  return(l)
}


results_list = list()
for(type in c("edge","path","node")){
  results_list[[type]] = list()
}

load_graph_vis_data = function(gsutil, outdir){
  
  files_to_read = c("gs://mawg-data/pass1b-06/merged/graphical_analysis_results_20220107.RData", # graphical results and pathway enrichments 
                    "gs://motrpac-data-freeze-pass/pass1b-06/v1.1/analysis/resources/motrpac-mappings-master_feature_to_gene.txt", # feature-to-gene map
                    "gs://mawg-data/pass1b-06/merged/kegg_reactome_pathway_hierarchies_20211221.RData") # pathway hierarchy lists 
  
  system(sprintf("mkdir -p %s", outdir))
  
  # check if files already exist
  target_files = sprintf("%s/%s", outdir, basename(files_to_read))
  files_to_download = files_to_read[!file.exists(target_files)]
  
  # download if necessary 
  if(length(files_to_download) > 0){
    system(sprintf("%s -m cp %s %s",
                   gsutil, 
                   paste(files_to_download, collapse=" "),
                   outdir))
  }
  
  # load RData
  load(sprintf("%s/%s", outdir, "graphical_analysis_results_20220107.RData"))
  load(sprintf("%s/%s", outdir, "kegg_reactome_pathway_hierarchies_20211221.RData"))
  
  # read in txt 
  feature_to_gene = fread(sprintf("%s/motrpac-mappings-master_feature_to_gene.txt", outdir), sep='\t', header=T)
  
  # process the enrichment analysis results
  fdr_enrichment_res = tree_analysis_enrichment_results[
    tree_analysis_enrichment_results$adj_p_value < 0.1 &
      tree_analysis_enrichment_results$intersection_size > 2,
  ]
  
  return(list(edge_sets = edge_sets,
              node_sets = node_sets, 
              feature_to_gene = as.data.frame(feature_to_gene),
              pathway_parents = pathway_parents,
              per_clust_enrichment_plots = per_clust_enrichment_plots, 
              repfdr_clusters = repfdr_clusters,
              repfdr_res = repfdr_res, 
              repfdr_clusters_pi = repfdr_clusters_pi,
              repfdr_clusters_str = repfdr_clusters_str, 
              sample_level_data = sample_level_data, 
              tree_analysis_cluster_df = tree_analysis_cluster_df,
              tree_analysis_enrichment_results = tree_analysis_enrichment_results,
              tree_analysis_selected_sets = tree_analysis_selected_sets,
              fdr_enrichment_res = fdr_enrichment_res, 
              zs_info = zs_info,
              zs_smoothed = zs_smoothed))
}

scratch='/tmp'
gsutil='/usr/local/anaconda3/bin/gsutil'
motrpac_mawg_path='/oak/stanford/groups/smontgom/nicolerg/src/MOTRPAC/motrpac-mawg'
feature_to_gene <- fread("~/data/smontgom/motrpac-mappings-master_feature_to_gene.txt", header = T, sep = "\t", fill = T)
load("~/data/smontgom/kegg_reactome_pathway_hierarchies_20211221.RData")
fdr_enrichment_res = tree_analysis_enrichment_results[
  tree_analysis_enrichment_results$adj_p_value < 0.1 &
    tree_analysis_enrichment_results$intersection_size > 2,
]
tissues = unique(cluster_res$tissue)


# load("~/data/smontgom/graphical_analysis_results_20211220.RData")
data = list(edge_sets = edge_sets,
            node_sets = node_sets, 
            feature_to_gene = as.data.frame(feature_to_gene),
            pathway_parents = pathway_parents,
            per_clust_enrichment_plots = per_clust_enrichment_plots, 
            repfdr_clusters = repfdr_clusters,
            repfdr_res = repfdr_res, 
            repfdr_clusters_pi = repfdr_clusters_pi,
            repfdr_clusters_str = repfdr_clusters_str, 
            sample_level_data = sample_level_data, 
            tree_analysis_cluster_df = tree_analysis_cluster_df,
            tree_analysis_enrichment_results = tree_analysis_enrichment_results,
            tree_analysis_selected_sets = tree_analysis_selected_sets,
            fdr_enrichment_res = fdr_enrichment_res, 
            zs_info = zs_info,
            zs_smoothed = zs_smoothed)
TISSUE <- "SKM-VL"
selected_clusters = data$tree_analysis_selected_sets[grepl(TISSUE, names(data$tree_analysis_selected_sets))]
clust <- names(selected_clusters)[[11]]

results_list = list()
for(type in c("edge","path","node")){
  results_list[[type]] = list()
}
for (clust in names(selected_clusters)){
  # what kind is it?
  # skip if all nodes are null 
  if(grepl("---", clust)){
    type = "edge"
    nodes = unname(unlist(strsplit(gsub(".*:","",clust), "---")))
    if(all(grepl("?w_F0_M0", nodes))){
      next()
    }
  }else if(grepl("->",clust)){
    type = "path"
    nodes = unname(unlist(strsplit(gsub(".*:","",clust), "->")))
    if(all(grepl("?w_F0_M0", nodes))){
      next()
    }
  }else{
    type = "node"
    if(grepl("?w_F0_M0", gsub(".*:","",clust))){
      next()
    }
  }
  
  # 1. tree
  tree = get_tree_plot_for_tissue(tissues = TISSUE,
                                  omes = NULL,
                                  node_sets = data$node_sets, 
                                  edge_sets = data$edge_sets,
                                  min_size = 10,
                                  highlight_subset = clust) 
  tree = tree + ggtitle(clust)
  
  # 2. grid of enrichments 
  enrich_res = data$tree_analysis_enrichment_results[data$tree_analysis_enrichment_results$cluster == clust,]
  enrich = plot_top_enrichments_for_tissue_set(
    enrich_res, 
    kegg_only=F,
    adj_pval_threshold=0.1, 
    max_num=20 # max number of pathways/rows in the plot
  )
  if(length(enrich$plot) > 0){
    enrich_grid = enrich$plot[[1]]
  }else{
    enrich_grid = NULL
  }
  
  # 3. sample-level trajectories 
  clust_data = data$tree_analysis_selected_sets[[clust]]
  sample_level = plot_group_mean_trajectories(
    clust_data,
    data$sample_level_data,
    gsutil = gsutil,
    motrpac_mawg_path = gitdir,
    scratch = scratch, 
    title = clust)
  
  # 4. network of enrichments (concatenated outside of this script)
  # visNetwork doesn't provide any way to print a plot within a loop
  # instead, save to HTML and concatenate HTML files after reports are knitted
  # to show multiple plots, we need to use a workaround that ONLY works when you KNIT TO HTML
  #system(sprintf("mkdir -p %s/visNetwork-html", outdir))
  #htmlout = sprintf("%s/visNetwork-html", outdir)
  tmpfile = sprintf("%s.html", clust) 
  tmpfile = gsub("->|:","__",tmpfile)
  vn = enrichment_network_vis(enrich_res,
                              data$feature_to_gene,
                              data$pathway_parents,
                              title = sprintf("%s, all datasets",clust),
                              out_html = tmpfile,
                              overwrite_html = T, # rerun if HTML is already generated
                              add_group_label_nodes = T,
                              return_html = T,
                              save_similarity_scores = sprintf("%s/%s", outdir, gsub("\\.html", "_similarity_scores\\.RDS", basename(tmpfile)))) 
  
  # save list of plots 
  results_list[[type]][[clust]] = list(tree = tree,
                                       top_enrich = enrich_grid,
                                       enrich_network = vn, 
                                       sample_level_traj = sample_level)
}

if(length(results_list$path)>0){
  cat("\n\n### Selected paths\n")
}
for(clust in names(results_list$path)){
  cat("\n\n#### ", clust,"\n")
  for(g in results_list$path[[clust]]){
    if(!is.null(g)){
      # print plot
      if(!is.character(g)){
        print(g)
        # create iframe
      }else if(endsWith(g, "html")){
        print(htmltools::tags$iframe(title = "My embedded document", 
                                     src = basename(g), 
                                     height = "960", width = "960"))
      }else{
        warning("I don't know what this is")
      }
    }
  }
  if(is.null(results_list$path[[clust]]$top_enrich)){
    cat(sprintf("No significant enrichments for %s\n", clust))
  }
  cat("\n\n---\n\n")
}

if(length(results_list$node)>0){
  cat("\n\n### Selected nodes\n")
}

if(length(results_list$edge)>0){
  cat("\n\n### Selected edges\n")
}
for(clust in names(results_list$edge)){
  cat("\n\n#### ", clust,"\n")
  for(g in results_list$edge[[clust]]){
    if(!is.null(g)){
      # print plot
      if(!is.character(g)){
        print(g)
        # create iframe
      }else if(endsWith(g, "html")){
        print(htmltools::tags$iframe(title = "My embedded document", 
                                     src = basename(g), 
                                     height = "960", width = "960"))
      }else{
        warning("I don't know what this is")
      }
    }
  }
  if(is.null(results_list$edge[[clust]]$top_enrich)){
    cat(sprintf("No significant enrichments for %s\n\n", clust))
  }
  cat("\n\n---\n\n")
}

#### testing function ####
















get_tree_plot_for_tissue<-function(
  tissues,
  omes = NULL,
  node_sets,
  edge_sets,
  min_size = 20,
  parallel_edges_by_ome = F,
  parallel_edges_by_tissue = F,
  edge_width_range = c(0,10),
  edge_alpha_range = c(0,1),
  color_nodes_by_states = T,
  max_trajectories = NULL,
  highlight_subset = NULL
){
  
  tissues = unique(tissues)
  omes = unique(omes)
  
  tissue_node_sets = limit_sets_by_regex(node_sets,tissues)
  tissue_node_sets = limit_sets_by_regex(tissue_node_sets,omes)
  tissue_edge_sets = limit_sets_by_regex(edge_sets,tissues)
  tissue_edge_sets = limit_sets_by_regex(tissue_edge_sets,omes)
  
  # Set things up to highlight a subset of the tree
  highlight_edge = F
  highlight_node = F
  highlight_path = F
  if(!is.null(highlight_subset)){
    if(grepl(":", highlight_subset)){
      # Check if it's this tissue 
      if(grepl(":", highlight_subset)){
        if(length(tissues)>1){
          warning("'highlight_subset' cannot handle a subset of tissues. Setting 'highlight_subset' to NULL Specify a single tissue or remove the tissue prefix from 'highlight_subset'.")
          highlight_subset = NULL
        }else{
          highlight_tissue = gsub(":.*","",highlight_subset)
          highlight_cluster = gsub(".*:","",highlight_subset)
          if(highlight_tissue != tissues){
            warning("'highlight_subset' tissue prefix does not match the tissue supplied. Setting 'highlight_subset' to NULL.")
            highlight_subset = NULL
          }
        }
      }else{
        highlight_cluster = highlight_subset
      }
    }
    if(!is.null(highlight_subset)){
      if(grepl("---", highlight_subset)){
        highlight_edge = T
      }else if(grepl("->",highlight_subset)){
        highlight_path = T
      }else{
        highlight_node = T
      }
    }
  }
  
  if(!is.null(highlight_subset) & any(parallel_edges_by_ome, parallel_edges_by_tissue, color_nodes_by_states)){
    warning("Setting 'parallel_edges_by_ome', 'parallel_edges_by_tissue', and 'color_nodes_by_states' to FALSE because 'highlight_subset' is TRUE.")
    parallel_edges_by_ome = F
    parallel_edges_by_tissue = F
    color_nodes_by_states = F
  }
  
  # Check the input omes and tissues sets
  # If either ome or tissues is NULL, take them from the data
  all_analytes = unique(unlist(node_sets))
  all_analytes_info = strsplit(all_analytes,split=";")
  data_tissue_set = unique(sapply(all_analytes_info,function(x)x[2]))
  data_ome_set = unique(sapply(all_analytes_info,function(x)x[1]))
  if(!all(tissues %in% data_tissue_set)){
    warning("Some input tissue names are not in the dataset and will be removed")
    tissues = intersect(tissues,data_tissue_set)
  }
  if(!all(omes %in% data_ome_set)){
    warning("Some input ome names are not in the dataset and will be removed")
    omes = intersect(omes,data_ome_set)
  }
  if(is.null(tissues)){
    tissues = data_tissue_set
  }
  if(is.null(omes)){
    omes = data_ome_set
  }
  
  # if we are not asked to get parallel edges and are asked to take 
  # the top trajectories then we can simply filter the edge sets and move on
  if(!is.null(max_trajectories) && !parallel_edges_by_tissue){
    tissue_edge_sets = filter_edge_sets_by_trajectories(tissue_edge_sets,max_trajectories)
  }
  
  # Create the initial matrices/data frame
  # with some additional info (depends on the user input), these will be used as input
  # for generating the graph
  node_info = cbind(
    names(tissue_node_sets),
    sapply(names(tissue_node_sets),function(x)strsplit(x,split="_")[[1]][1]),
    sapply(tissue_node_sets,length)
  )
  colnames(node_info) = c("node","week","size")
  edge_info = cbind(
    sapply(names(tissue_edge_sets),function(x)strsplit(x,split="---")[[1]][1]),
    sapply(names(tissue_edge_sets),function(x)strsplit(x,split="---")[[1]][2]),
    sapply(tissue_edge_sets,length)
  )
  colnames(edge_info) = c("prev","next","size")
  
  # If we are asked to add parallel edges by omes then we need to 
  # compute ome edge sizes and add them
  if(parallel_edges_by_ome){
    for(o in omes){
      v = c()
      for(e in rownames(edge_info)){
        v[e] = sum(grepl(paste0(o,";"),tissue_edge_sets[[e]]))
      }
      edge_info = cbind(edge_info,v)
      colnames(edge_info)[ncol(edge_info)] = o
    }
    # always use the same assay colors
    edge_colors = define_assay_cols() # defined in cluster_viz_fx.R
    # require(RColorBrewer)
    # o = c("METAB","TRNSCRPT","PROT","ACETYL","PHOSPHO","UBIQ","ATAC","METHYL","IMMUNO")
    # edge_colors = brewer.pal(length(o), 'Set1')
    # names(edge_colors) = o
  }
  # if parallel_edges_by_tissue = TRUE then:
  # if we do not need to filter the top trajectories then we can simply count
  # the number of features per edge for a tissue using regular expressions
  if(is.null(max_trajectories) && parallel_edges_by_tissue){
    for(tissue in tissues){
      v = c()
      for(e in rownames(edge_info)){
        v[e] = sum(grepl(paste0(tissue,";"),tissue_edge_sets[[e]]))
      }
      edge_info = cbind(edge_info,v)
      colnames(edge_info)[ncol(edge_info)] = tissue
    }
    # always use the same assay colors
    edge_colors = MotrpacBicQC::tissue_cols
  }
  # if parallel_edges_by_tissue = TRUE then:
  # if we are asked to filter by the top trajectories then we need to 
  # recompute the edge sizes per tissue by limiting the tissue edge sets
  if(!is.null(max_trajectories) && parallel_edges_by_tissue){
    for(tissue in tissues){
      curr_tissue_edge_sets = limit_sets_by_regex(tissue_edge_sets,tissue)
      curr_tissue_edge_sets = filter_edge_sets_by_trajectories(
        curr_tissue_edge_sets,max_trajectories)
      v = c()
      for(e in rownames(edge_info)){
        v[e] = length(curr_tissue_edge_sets[[e]])
      }
      edge_info = cbind(edge_info,v)
      colnames(edge_info)[ncol(edge_info)] = tissue
    }
    # always use the same assay colors
    edge_colors = MotrpacBicQC::tissue_cols
  }
  
  # transform the initial node/edge info matrices to data frames
  # filter edge sizes by the specified user input
  d = as.data.frame(edge_info)
  for(j in 3:ncol(d)){
    d[[j]] = as.numeric(d[[j]])
    d[d[,j]<min_size,j] = 0 # filter by min edge size
  }
  d_nodes = data.frame(node_info,check.names = F)
  d_nodes[[3]] = as.numeric(d_nodes[[3]])
  d_nodes$inds = 0:(nrow(d_nodes)-1)
  
  # If we need to plot parallel edges then we need to transform the edge-based
  # data frame into a long format (to make it work with ggraph)
  if(parallel_edges_by_ome || parallel_edges_by_tissue){
    added_cols = names(d)[-c(1:3)]
    if(!parallel_edges_by_ome){added_cols = setdiff(added_cols,omes)}
    if(!parallel_edges_by_tissue){added_cols = setdiff(added_cols,tissues)}
    newd = d[,1:2]
    newd$size = 0
    newd$type = NA
    for(o in added_cols){
      currd = d[,c(1:2,which(colnames(d)==o))]
      names(currd)[3] = "size"
      currd$type = o
      newd = rbind(newd,currd)
    }
    newd = newd[!is.na(newd$type),]
    d_g = igraph::graph_from_data_frame(newd)
    if(all(E(d_g)$size == 0)){
      message("No non-0 edges. Skipping.")
      return()
    }
    E(d_g)$edge_size = edge_attr(d_g,"size")
    d_g_auto_layout <- create_layout(d_g, layout = 'auto')
  }else{
    d_g = igraph::graph_from_data_frame(d)
    if(all(E(d_g)$size == 0)){
      message("No non-0 edges. Skipping.")
      return()
    }
    E(d_g)$edge_size = edge_attr(d_g,"size")
    d_g_auto_layout <- create_layout(d_g, layout = 'auto')
  }
  
  ################
  # From this point we start generating ggraph info
  # (e.g., layout, labels etc.)
  d_g_ordered_nodes = c("F-1_M1","F1_M-1","F0_M-1","F-1_M-1","F-1_M0",
                        "F0_M0","F0_M1","F1_M1","F1_M0")
  d_g_ordered_cols = c("lightgreen","lightgreen","blue","blue","blue",
                       "gray","red3","red3","red3")
  #d_g_ordered_cols_alt = c("black","black","black","black","black",
  #                     "gray","black","black","black")
  d_g_ordered_cols_alt = c("gray4","gray4","gray4","gray4","gray4",
                           "gray","gray4","gray4","gray4")
  d_g_ordered_shapes = c(18,18,19,19,19,19,19,19,19)
  layer_plot_names = c(
    "F down, M up","F up, M down",
    "M only down","Both down","F only down",
    "No\nchange",
    "M only up","Both up","F only up"
  )
  
  names(d_g_ordered_cols) = d_g_ordered_nodes
  names(d_g_ordered_shapes) = d_g_ordered_nodes
  names(d_g_ordered_cols_alt) = d_g_ordered_nodes
  l_x_lim = c(min(d_g_auto_layout$x),max(d_g_auto_layout$x))
  l_y_lim = c(min(d_g_auto_layout$y),max(d_g_auto_layout$y))
  d_g_our_layout  = copy(d_g_auto_layout)
  # set 5 horiz layers over time
  xjump = (l_x_lim[2]-l_x_lim[1]) / 4
  # set 9 vertical layers over ordered_nodes
  yjump = (l_y_lim[2]-l_y_lim[1]) / 8
  d_g_our_layout[d_g_our_layout$name=="0w","x"] = l_x_lim[1]
  d_g_our_layout[d_g_our_layout$name=="0w","y"] = mean(l_y_lim)
  curr_weeks = c("1w","2w","4w","8w")
  for(j in 1:length(curr_weeks)){
    w = curr_weeks[j]
    d_g_our_layout[grepl(w,d_g_our_layout$name),"x"] = l_x_lim[1] + j*xjump
  }
  for(j in 1:length(d_g_ordered_nodes)){
    n = d_g_ordered_nodes[j]
    d_g_our_layout[grepl(n,d_g_our_layout$name),"y"] = l_y_lim[1] + (j-1)*yjump
  }
  # make sure that the 0w node is in the same line as of the no response
  d_g_our_layout[1,"y"] = d_g_our_layout[grepl("F0_M0",d_g_our_layout$name),"y"][1]
  # set node sizes and other features
  V(d_g)$setsize = d_nodes[V(d_g)$name,"size"]
  V(d_g)$setsize[V(d_g)$name == "0w"] = median(V(d_g)$setsize)
  V(d_g)$label = sapply(V(d_g)$name,
                        function(x){a=strsplit(x,split="w_")[[1]];a[length(a)]})
  V(d_g)$label = gsub("_","\n",V(d_g)$label)
  V(d_g)$col = "gray"
  V(d_g)$alt_col = "gray"
  V(d_g)$shape = 15
  for(j in 1:length(d_g_ordered_nodes)){
    n = d_g_ordered_nodes[j]
    V(d_g)$col[grepl(n,d_g_our_layout$name)] = d_g_ordered_cols[n]
    V(d_g)$alt_col[grepl(n,d_g_our_layout$name)] = d_g_ordered_cols_alt[n]
    V(d_g)$shape[grepl(n,d_g_our_layout$name)] = d_g_ordered_shapes[n]
  }
  # for(ome in colnames(tissue_ome_data)){
  #   d_g = set_vertex_attr(d_g,ome,V(d_g),tissue_ome_data[V(d_g)$name,ome])
  # }
  grid_group_annotation_y = c()
  for(j in 1:length(d_g_ordered_nodes)){
    n = d_g_ordered_nodes[j]
    grid_group_annotation_y[n] = l_y_lim[1] + (j-1)*yjump
  }
  weeks_annotation_x = c()
  for(j in 1:length(curr_weeks)){
    w = curr_weeks[j]
    weeks_annotation_x[w] = l_x_lim[1] + j*xjump
  }
  
  # Set colors to highlight edge, node, path
  # Red to highlight; gray for everything else 
  if(highlight_edge){
    V(d_g)$col = "gray"
    V(d_g)$alt_col = "gray"
    # Make edge red
    ecols = rep("gray", nrow(d))
    ecols[which(rownames(d)==highlight_cluster)] = "red"
    E(d_g)$col = ecols
  }else if(highlight_node){
    E(d_g)$col = "gray"
    vcols = rep("gray", length(names(V(d_g))))
    vcols[which(names(V(d_g))==highlight_cluster)] = "red"
    V(d_g)$col = vcols
    V(d_g)$alt_col = vcols
  }else if(highlight_path){
    # extract nodes
    curr_nodes = c("0w",unname(unlist(strsplit(highlight_cluster, "->"))))
    # extract edges
    edge_vector = c()
    for(e in 2:length(curr_nodes)){
      edge_vector = c(edge_vector, sprintf("%s---%s", curr_nodes[e-1], curr_nodes[e]))
    }
    # set colors for edges
    ecols = rep("gray", nrow(d))
    ecols[rownames(d)%in%edge_vector] = "red"
    E(d_g)$col = ecols
    # set colors for nodes
    vcols = rep("gray", length(names(V(d_g))))
    vcols[names(V(d_g))%in%curr_nodes] = "red"
    V(d_g)$col = vcols
    V(d_g)$alt_col = vcols
  }
  
  if(parallel_edges_by_ome || parallel_edges_by_tissue){
    p = ggraph(d_g,layout = d_g_our_layout,) +  
      geom_edge_fan0(
        aes(width = E(d_g)$edge_size,alpha=E(d_g)$edge_size,colour = E(d_g)$type),
        end_cap = circle(3, 'mm'),curvature = 0.1) +
      scale_edge_color_discrete(name="")
  }else if(!is.null(highlight_subset)){
    #CHANGED A BIT HERE
    p = ggraph(d_g,layout = d_g_our_layout,) +  
      geom_edge_arc(aes(width = E(d_g)$edge_size, alpha=E(d_g)$edge_size, 
                        color = rep("gray", length(E(d_g)$col))), #make all the initial edges grey
                    end_cap = circle(3, 'mm'),curvature = 0.1) +  
      scale_edge_colour_identity(guide = "none")
  }else{
    p = ggraph(d_g,layout = d_g_our_layout,) +  
      #geom_edge_link(aes(width = E(d_g)$edge_size,alpha=E(d_g)$edge_size)) + 
      geom_edge_arc(aes(width = E(d_g)$edge_size,alpha=E(d_g)$edge_size),
                    end_cap = circle(3, 'mm'),curvature = 0.1) +
      scale_edge_color_discrete(name="")
  }
  
  p = p +
    scale_edge_width(range=edge_width_range,name="Intersect size") + 
    scale_edge_alpha(range=edge_alpha_range,name="Intersect size") 
  
  if(color_nodes_by_states){
    p = p + geom_node_point(aes(size = V(d_g)$setsize),alpha=1,
                            color=V(d_g)$col,shape=V(d_g)$shape)
  }else{
    #STARTED CHANGING STUFF HERE
    rightproper_node_sizes <- rep(NA, length(V(d_g)$setsize)) #don't want to plot other nodes cos then you double opacity
    rightproper_node_sizes[names(V(d_g))%in%curr_nodes] <- 50 #set the size of the path-set here
    
    rightproper_edge_widths <- rep(NA, length(E(d_g)$edge_size)) #same thing with edges widths
    rightproper_edge_widths[rownames(d)%in%edge_vector] <- 50 #set the size of the path-set here
    rightproper_edge_alphas <- rep(NA, length(E(d_g)$edge_size)) #same thing with alphas
    rightproper_edge_alphas[rownames(d)%in%edge_vector] <- 100 #set the size of the path-set here, or just make them max opacity
    p = p + 
      geom_node_point(aes(size = V(d_g)$setsize),alpha=1,
                      color=rep("gray", length(V(d_g)$alt_col)),shape=V(d_g)$shape) + 
      geom_node_point(aes(size = rightproper_node_sizes),alpha=1, #overlays white nodes on grey ones, don't need this if fully opaque
                      color="white",shape=V(d_g)$shape) + 
      geom_node_point(aes(size = rightproper_node_sizes),alpha=1, #overlays red circles for nodes
                      color=V(d_g)$alt_col,shape=V(d_g)$shape) +  
      geom_edge_arc(aes(width = rightproper_edge_widths, alpha=max(E(d_g)$edge_size), color = "white"), #overlays white edges on grey ones
                    end_cap = circle(3, 'mm'),curvature = 0.1) +
      geom_edge_arc(aes(width = rightproper_edge_widths, alpha=max(E(d_g)$edge_size), color = E(d_g)$col), #overlays red edges on white ones
                    end_cap = circle(3, 'mm'),curvature = 0.1) +
      scale_edge_colour_identity(guide = "none")
  }
  p = p +
    #geom_node_text(aes(label = V(d_g)$label), repel=F) + 
    scale_size(range = c(2,20),name="Number of analytes") +
    theme(panel.background = element_rect(fill="white"),
          legend.key.size = unit(0.4, 'cm'))
  
  names(grid_group_annotation_y) = layer_plot_names
  for(g in names(grid_group_annotation_y)){
    p = p + annotate(geom="text",x=l_x_lim[1] + 0.2, y=grid_group_annotation_y[g],
                     label=g,fontface = "bold")
  }
  names(weeks_annotation_x) = paste("week",c(1,2,4,8))
  for(w in names(weeks_annotation_x)){
    p = p + annotate(geom="text",y=l_y_lim[2]+0.5, x=weeks_annotation_x[w],
                     label=w,fontface = "bold")
  }
  #plot(p)
  
  # change colors and legend edge size
  if(parallel_edges_by_ome | parallel_edges_by_tissue){
    p = p +
      scale_edge_color_manual(values=edge_colors, name="") +
      guides(edge_color=guide_legend(override.aes = list(edge_width=3)))
  }
  
  return(p)
}


tree = get_tree_plot_for_tissue(tissues = TISSUE,
                                omes = NULL,
                                node_sets = data$node_sets, 
                                edge_sets = data$edge_sets,
                                min_size = 10,
                                highlight_subset = clust) 
tree = tree + ggtitle(clust)
tree

