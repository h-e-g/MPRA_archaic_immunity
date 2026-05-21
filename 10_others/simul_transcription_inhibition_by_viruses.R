MPRA_DIR <- "/pasteur/helix/projects/evo_immuno_pop/MPRA"
ANALYSIS_DIR <- "MPRA_count_exp6_analysis3"

source(sprintf("%s/scripts/misc_plots.R", MPRA_DIR))

per_cell_integ=round(rnorm(1000,200,10))

## 4 transcripts per BC on average
per_cell_UMI=round(rnorm(1000,8000,10))

#barcodes numbers from 1:1000
perBC_logActivity=rnorm(1000,0,1)
perBC_Activity=2^perBC_logActivity

perCell_BC=lapply(per_cell_integ, function(x) sample(1:1000,x,replace=FALSE))
perCell_counts=lapply(seq_along(perCell_BC), function(i){
  BCs=perCell_BC[[i]]
  UMI_perBC=sample(BCs,prob=perBC_Activity[BCs]/sum(perBC_Activity[BCs]),size=per_cell_UMI[i],replace=TRUE)
  tab=table(UMI_perBC)
  data.table(cell=i,count=as.numeric(tab),BC=as.numeric(names(tab)))
  })
perCell_counts = rbindlist(perCell_counts)
perCell_counts[,length(unique(BC)),by=cell]
perCell_counts[,sum(count),by=cell]
UMIs_normal=perCell_counts[,.(DNA=length(unique(cell)),RNA=sum(count)),by=BC]



## half of cells infected with reduced transcription
per_cell_UMI_infected=c(round(rnorm(200,8000,10)),pmax(round(rnorm(800,20,3)),0))

perCell_counts_infected=lapply(seq_along(perCell_BC), function(i){
  BCs=perCell_BC[[i]]
  UMI_perBC=sample(BCs,prob=perBC_Activity[BCs]/sum(perBC_Activity[BCs]),size=per_cell_UMI_infected[i],replace=TRUE)
  tab=table(UMI_perBC)
  data.table(cell=i,count=as.numeric(tab),BC=as.numeric(names(tab)))
  })
perCell_counts_infected = rbindlist(perCell_counts_infected)
UMIs_infected=perCell_counts_infected[,.(DNA=length(unique(cell)),RNA=sum(count)),by=BC]


UMIs_all=merge(UMIs_normal,UMIs_infected,by='BC',suffix=c('_normal','_infected'))
UMIs_all[,log2_activity_normal:=log2((RNA_normal/DNA_normal)/(sum(RNA_normal)/sum(DNA_normal)))]
UMIs_all[,log2_activity_infected:=log2((RNA_infected/DNA_infected)/(sum(RNA_infected)/sum(DNA_infected)))]
UMIs_all[,log2FC:=log2_activity_infected-log2_activity_normal]

library(ggplot2)
library(ggrastr)
p <- ggplot(UMIs_all,aes(x=log2_activity_normal, y=log2_activity_infected)) 
p <- p + geom_vline(xintercept=0,col='lightgrey',linetype=2) + geom_hline(yintercept=0,col='lightgrey',linetype=2)
p <- p + geom_abline(intercept=0,slope=1,col='darkgrey')
p <- p + rasterize(geom_point(size=.4),dpi=400)
p <- p + xlab('log2 CRS activity (basal)') + ylab('log2 CRS activity (infected)')
p <- p + theme_plot()

pdf(sprintf('%s/figures/%s/Simulation_shutdown.pdf',MPRA_DIR, ANALYSIS_DIR),width=3,height=3)
print(p)
dev.off()

# To show that a global shutdown of transcription in a subset of cells is enough to mimic the patter we observed in viral conditions
# a performed simulations where N=1000 BCs are integrated across 1000 cells, with an average ofr 200 BC per cell
# each barcode has log-normally distributed activity (ie. normalized RNA/DNA ratio) 
# each cell transcribed a fixed amount of UMIs either normaly distributed aroudn 8000 (normal or bysrander cells), or 20 (infected cells)
# in the infected condition only a fixed percenntage of cells are infected.


# what to show in supplementary note.
# show : BC activity decreasing order %o% UMI per cells
# cells ordered by UMIs

# RNA per cell normal  # RNA per cell infected   # DNA per cell,
# activity ratios normal, activity ratios infected # normal VS infected

