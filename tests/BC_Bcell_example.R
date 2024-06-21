# library(SGR)

library(devtools)
load_all()

BC_Bcells_u <- readRDS("data/BC_B_cells_u.rds")
BC_loc <- readRDS("data/BCIDC_location.RDS")
BC_loc[,1] <- (BC_loc[,1] - mean(BC_loc[,1])) / sd(BC_loc[,1])
BC_loc[,2] <- (BC_loc[,2] - mean(BC_loc[,2])) / sd(BC_loc[,2])
BC_loc <- as.matrix(BC_loc)

n <- dim(BC_Bcells_u)[1]
p <- dim(BC_Bcells_u)[2]

#BC_Bcell_est <- SGR_undirected(BC_Bcells_u, BC_loc, p, 100, 100, 1)
#list.save(BC_Bcell_est,"BC_Bcell_est.rds")

abc <- SGR_undirected(BC_Bcells_u[,1:4], BC_loc, 4, 5, 5, 1)
dir.create('./example_output')
list.save(abc,"example_output/BC_Bcell_est_example.rds")
