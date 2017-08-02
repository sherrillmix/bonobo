all: getData.md makeOtus.md blastData.md parseBlast.md plotPcoa.md plotHeatmap.md 16s/runQiime.md 16s/plotShannon.md 16s/plotHeatmap.md 16s/plotPcoa.md 16s/plotBetaDiversity.md sensitivity/sensitivity.md README.md README.html

getData.md: getData.Rmd
	R -e 'knitr::knit("getData.Rmd")'

makeOtus.md: makeOtus.Rmd getData.md
	R -e 'knitr::knit("makeOtus.Rmd")'

blastData.md: blastData.Rmd makeOtus.md
	R -e 'knitr::knit("blastData.Rmd")'

parseBlast.md: parseBlast.Rmd blastData.md
	R -e 'knitr::knit("parseBlast.Rmd")'

plotPcoa.md: plotPcoa.Rmd parseBlast.md
	R -e 'knitr::knit("plotPcoa.Rmd")'

plotHeatmap.md: plotHeatmap.Rmd parseBlast.md
	R -e 'knitr::knit("plotHeatmap.Rmd")'
	
16s/runQiime.md: 16s/runQiime.Rmd getData.md
	cd 16s && R -e 'knitr::knit("runQiime.Rmd")'

16s/plotShannon.md: 16s/plotShannon.Rmd 16s/runQiime.md
	cd 16s && R -e 'knitr::knit("plotShannon.Rmd")'

16s/plotHeatmap.md: 16s/plotHeatmap.Rmd 16s/runQiime.md
	cd 16s && R -e 'knitr::knit("plotHeatmap.Rmd")'

16s/plotPcoa.md: 16s/plotPcoa.Rmd 16s/runQiime.md
	cd 16s && R -e 'knitr::knit("plotPcoa.Rmd")'

16s/plotBetaDiversity.md: 16s/plotBetaDiversity.Rmd 16s/runQiime.md 16s/plotPcoa.md
	cd 16s && R -e 'knitr::knit("plotBetaDiversity.Rmd")'

sensitivity/sensitivity.md: sensitivity/sensitivity.Rmd
	cd sensitivity && R -e 'knitr::knit("sensitivity.Rmd")'

README.md: README.template getData.md makeOtus.md blastData.md parseBlast.md plotPcoa.md plotHeatmap.md 16s/runQiime.md 16s/plotShannon.md 16s/plotHeatmap.md 16s/plotPcoa.md 16s/plotBetaDiversity.md sensitivity/sensitivity.md
cat README.template getData.md makeOtus.md blastData.md parseBlast.md plotPcoa.md plotHeatmap.md <( sed 's@(figure/@(16s/figure/@' 16s/runQiime.md 16s/plotShannon.md 16s/plotHeatmap.md 16s/plotPcoa.md 16s/plotBetaDiversity.md) <( sed 's@(figure/@(sensitivity/figure/@' sensitivity/sensitivity.md) > README.md

README.html: README.md
	pandoc --standalone --smart --self-contained  --css=github-pandoc.css --toc README.md -o README.html

