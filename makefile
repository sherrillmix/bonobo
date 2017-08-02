all: getData.md makeOtus.md blastData.md parseBlast.md makeClusters.md plotHeatmap.md 16s/runQiime.md 16s/plotShannon.md 16s/plotHeatmap.md 16s/plotPcoa.md 16s/plotBetaDiversity.md

getData.md: getData.Rmd
	R -e 'knitr::knit("getData.Rmd")'

makeOtus.md: makeOtus.Rmd getData.Rmd
	R -e 'knitr::knit("makeOtus.Rmd")'

blastData.md: blastData.Rmd makeOtus.md
	R -e 'knitr::knit("blastData.Rmd")'

parseBlast.md: parseBlast.Rmd blastData.md
	R -e 'knitr::knit("parseBlast.Rmd")'

makeClusters.md: makeClusters.Rmd parseBlast.md
	R -e 'knitr::knit("makeClusters.Rmd")'

plotHeatmap.md: plotHeatmap.Rmd parseBlast.md
	R -e 'knitr::knit("plotHeatmap.Rmd")'
	
16s/runQiime.md: 16s/runQiime.Rmd getData.Rmd
	cd 16s && R -e 'knitr::knit("runQiime.Rmd")'

16s/plotShannon.md: 16s/plotShannon.Rmd 16s/runQiime.md
	cd 16s && R -e 'knitr::knit("plotShannon.Rmd")'

16s/plotHeatmap.md: 16s/plotHeatmap.Rmd 16s/runQiime.md
	cd 16s && R -e 'knitr::knit("plotHeatmap.Rmd")'

16s/plotPcoa.md: 16s/plotPcoa.Rmd 16s/runQiime.md
	cd 16s && R -e 'knitr::knit("plotPcoa.Rmd")'

16s/plotBetaDiversity.md: 16s/plotBetaDiversity.Rmd 16s/runQiime.md 16s/plotPcoa.md
	cd 16s && R -e 'knitr::knit("plotBetaDiversity.Rmd")'
