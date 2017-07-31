all: getData.md makeOtus.md blastData.md 16s/runQiime.md parseBlast.md

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


	
16s/runQiime.md: 16s/runQiime.Rmd getData.Rmd
	cd 16s && R -e 'knitr::knit("runQiime.Rmd")'

16s/shannonOtus.md: 16s/shannonOtus.Rmd 16s/runQiime.md
	cd 16s && R -e 'knitr::knit("shannonOtus.Rmd")'
