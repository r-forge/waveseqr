checkFolderHierarchy <- function(indir){
#
#   Checks the presence of directories 'log', 'bed', 'graph', 'padgraph',
#   'thresdist', 'peaks' and creates them if absent
#

#   Apratim Mitra 2011-2012

 folders <- c("log","bed","graph","padgraph","thresdist","peaks")
 for(i in 1:length(folders)){    
    dir.create(file.path(indir,folders[i]),showWarnings=FALSE)
 }
}