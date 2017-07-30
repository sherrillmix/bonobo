all: getData.md makeOtus.md 16s/runQiime.md

getData.md: getData.Rmd
	R -e 'knitr::knit("getData.Rmd")'

makeOtus.md: makeOtus.Rmd getData.Rmd
	R -e 'knitr::knit("makeOtus.Rmd")'
	
16s/runQiime.md: 16s/runQiime.Rmd getData.Rmd functions.R
	R -e 'knitr::knit("16s/runQiime.Rmd")'

