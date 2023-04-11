library(SpATS)
library(data.table)

setwd("~/Documents/git/Roza_2021/figs/")

a1 <- read.csv("~/Documents/git/Roza_2021/raw_data/Roza 2021_yield.csv")
a1 <- a1[,-c(1,2)]
lev2 <- c("gen","rep")
a1[,lev2] <- lapply(a1[,lev2], factor)
nlevels(a1$gen)
nlevels(a1$col)

a1$R <- as.factor(a1$row)
a1$C <- as.factor(a1$col)
nlevels(a1$C)
nlevels(a1$R)
head(a1)
colnames(a1)

lev3 <- colnames(a1)[5:14]

Y1 <- list()
H2 <- list()
for (i in 1:length(lev3)) {
  m1 <- SpATS(response = lev3[i],
              spatial = ~PSANOVA(col, row, nseg = c(nlevels(a1$C),nlevels(a1$R)), degree = c(3,3), nest.div = 2),
              fixed = ~rep,
              random = ~ rep:C + rep:R,
              genotype.as.random = T,
              genotype = "gen",
              data = a1,
              control = list(tolerance = 1e-03))
  
  h2 <- getHeritability(m1)
  h2 <- as.data.frame(h2)
  H2[[length(H2)+1]] <- h2
  
  # pred4 <- predict.SpATS(m1, which = "gen", predFixed = "marginal")
  # pred4 <- pred4[,c(1,7,8)]
  # pred4$weight <- (1/pred4$standard.errors)^2
  # Y1[[length(Y1)+1]] <- pred4
  # 
}
names(Y1) <- lev3
names(H2) <- lev3


Y2 <-rbindlist(Y1, use.names=TRUE, fill=TRUE, idcol="env")
H2 <-rbindlist(H2, use.names=TRUE, fill=TRUE, idcol="env")

write.csv(H2, "H2.csv", row.names = F, quote = F)

plot(m1)
m1

head(Y2)
Y2$env <- as.factor(Y2$env)
levels(Y2$env)
lev6 <- c("sep_21_Y","may_22_Y","jul_22_Y","aug_22_Y","sep_22_Y")

Y4 <- Y2 %>% dplyr::filter(env %in% lev6)
Y4 <- droplevels(Y4)
levels(Y4$env)
head(Y4)
Y4 <- Y4 %>% separate(1, c("month", "year","trait"), sep = "_", remove = F, convert = FALSE, extra = "merge")
lev5 <- c("month", "year")
Y4 <- as.data.frame(Y4)
Y4 <- Y4[,-4]
Y4[,lev5] <- lapply(Y4[,lev5], factor)
Y4$month <- fct_relevel(Y4$month, c("may", "jul", "aug", "sep"))
str(Y4)

levels(Y4$month)
levels(Y4$year)
levels(Y4$gen)

data <- Y4
data <- data[order(data$gen, data$year), ]
data1 <- na.omit(data)
head(data1)
str(data1)

FA_1 <- asreml::asreml(fixed = predicted.values ~ 1, 
                       random = ~ fa(year, 1):ar1(month):id(gen) + env,
                       data = data1, na.action = list(x = "include", y = "include"), 
                       weights = weight, family = asreml::asr_gaussian(dispersion = 1))


FA_1 <- update.asreml(FA_1)
summary(FA_1)$varcomp
current.asrt <- as.asrtests(FA_1, NULL, NULL)
current.asrt <- rmboundary.asrtests(current.asrt)

# asreml.options(workspace="1024mb")
# this part collapse the program
diffs <- predictPlus(classify = "gen:env",
                     asreml.obj = FA_1,
                     wald.tab = current.asrt$wald.tab,
                     present = c("month","gen","year","env"))

ST1 <-diffs[[1]]


diffs2 <- predictPlus(classify = "gen:month", 
                     asreml.obj = FA_1, 
                     wald.tab = current.asrt$wald.tab, 
                     present = c("month","gen","year","env"))
ST2 <-diffs2[[1]]

diffs3 <- predictPlus(classify = "gen:year", 
                     asreml.obj = FA_1, 
                     wald.tab = current.asrt$wald.tab, 
                     present = c("month","gen","year","env"))
ST3 <-diffs3[[1]]

diffs4 <- predictPlus(classify = "gen", 
                     asreml.obj = FA_1, 
                     wald.tab = current.asrt$wald.tab, 
                     present = c("month","gen","year","env"))
ST4 <-diffs4[[1]]

head(ST1)
head(ST2)
head(ST3)
head(ST4)


ST2$month <- paste0("ST2_", ST2$month)
ST3$year <- paste0("ST3_", ST3$year)

colnames(ST4)[2] <- "ST4_Yi"
ST4 <- ST4 %>% select(1:2)

ST2 <- ST2 %>% select(1:3) %>% spread(key = month, value = predicted.value) 
ST3 <- ST3 %>% select(1:3) %>% spread(key = year, value = predicted.value)
ST2 <- ST2[,c(1,4,3,2,5)]

Y2.1 <- Y2
Y2.1$env <- paste0("ST0_", Y2.1$env)
Y2.1$env <- as.factor(Y2.1$env)
levels(Y2.1$env)
ST0 <- Y2.1 %>% select(1:3) %>% spread(key = env, value = predicted.values)


head(ST1)
ST1$env <- paste0("ST1_", ST1$env)
ST1 <- ST1 %>% select(1:3) %>% spread(key = env, value = predicted.value)

BLUP5 <- inner_join(ST0, ST1, by = "gen") %>% inner_join(., ST2, by = "gen") %>% inner_join(., ST3, by = "gen")  %>% inner_join(., ST4, by = "gen")
colnames(BLUP5)

BLUP6 <- BLUP5 %>% dplyr::select(c("gen",
                                   "ST0_sep_21_Y","ST0_may_22_Y","ST0_jul_22_Y","ST0_aug_22_Y","ST0_sep_22_Y",
                                   "ST1_sep_21_Y","ST1_may_22_Y","ST1_jul_22_Y","ST1_aug_22_Y","ST1_sep_22_Y",
                                   "ST2_may","ST2_jul","ST2_aug","ST2_sep","ST3_21","ST3_22","ST4_Yi"))

colnames(BLUP6)

write.csv(BLUP6, "~/Documents/git/Roza_2021/raw_data/BLUEs_Yi_Roza2021.csv", quote = F, row.names = F)

rm(FA_1)
rm(current.asrt)

# save.image("~/Documents/git/big_files/FA1.RData")
# load("~/Documents/git/big_files/FA1.RData")
