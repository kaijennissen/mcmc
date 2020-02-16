devtools::install_github("csgillespie/efficient",
						 args = "--with-keep.source")

library(tidyverse)
library(efficient)
library(foreach)
library(doParallel)

cl <- parallel::makeForkCluster(4)
doParallel::registerDoParallel(cl)
monopoly_resu = 

tic()
monopoly_resu = foreach(i = 1:1000) %dopar% {
	monopoly_resu = simulate_monopoly(10000)
}
toc()


tbl <- reduce(monopoly_resu , cbind)
tbl <- t(tbl)
colnames(tbl) <- paste0("field_", c(1:40))
rownames(tbl) <- paste0("run_", c(1:1000))
tbl <- as_tibble(tbl) %>% 
	gather(field, freq) %>%
	separate(field,c( "name", "field")) %>% select(-name) %>% 
	mutate(field=as.factor(field))

ggplot(tbl, aes(x=field, y=freq))+
	geom_boxplot()

tibble(field=as.factor(1:40), mean = map_dbl(tbl, mean), 
	   q25 = map_dbl(tbl, ~quantile(.x, .25)),
	   q75 = map_dbl(tbl, ~quantile(.x, .75))
	   )

%>% 
ggplot()+
	geom_ribbon(aes(ymin=q25, ymax=q75, x=field), fill = "lightgray")+
	geom_point(aes(x = field, y = mean))





