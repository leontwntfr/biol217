#--------------------------

data("midwest")

## popdensity in state (boxplot)

popdens <- ggplot(midwest, aes(state, popdensity, fill = state)) + geom_boxplot(coef = 6)

popdens + ylim(min = 0, max = 5000) +
  scale_fill_hue(name = "State", labels = c("Illinois", "Indiana", 'Michigan', 'Ohio', 'Wisconsin')) + 
  labs(y = 'Population density', x = 'State', title = 'Comparison of Population Densities Among Midwestern States')

ggsave('popdens.png', height = 6, width = 8, unit = 'in', dpi = 300)



## pop in state (boxplot)

pop <- ggplot(midwest, aes(state, poptotal, fill = state)) + geom_boxplot(coef = 6)

pop + ylim(min = 0, max = 100000) +
  scale_fill_hue(name = "State", labels = c("Illinois", "Indiana", 'Michigan', 'Ohio', 'Wisconsin')) + 
  labs(y = 'Population', x = 'State', title = 'Comparison of Total Populations Among Midwestern States')

ggsave('pop.png', height = 6, width = 8, unit = 'in', dpi = 300)


## county count

county <- ggplot(midwest, aes(state)) + geom_bar() + 
  scale_x_discrete(labels= c("Illinois", "Indiana", 'Michigan', 'Ohio', 'Wisconsin'))

ggsave('counties.png', height = 6, width = 8, unit = 'in', dpi = 300)

## poverty - college

pov_college <- ggplot(midwest, aes(percollege, percbelowpoverty, col = state)) + geom_point() + 
  scale_color_discrete(name = "State", labels = c("Illinois", "Indiana", 'Michigan', 'Ohio', 'Wisconsin')) +
  geom_smooth(se = FALSE) +
  labs(y = 'Percentage of population below poverty line', x = 'Percentage of population with college education', title = 'Below Poverty Line - College Education Plot for Midwestern States')

ggsave('pover_college.png', height = 6, width = 8, unit = 'in', dpi = 300)




midwest %>% select(state, poptotal, popwhite, popblack, popamerindian, popasian, popother) -> subsetmidwest

illi <- subsetmidwest[1:102,]

popsum <- data.frame(sum(illi$popwhite), sum(illi$popblack), sum(illi$popamerindian), sum(illi$popasian), sum(illi$popother))




#----------------------------------

data('chickwts')

ggplot(chickwts, aes(feed, weight)) + geom_boxplot(coef = 6) + 
  labs(y = 'Weight of chicks [g]', x = 'Feed type', title = 'Comparison of Weights of Chicks Given Different Feed After Six Weeks') + 
  scale_x_discrete(labels= c("Casein", "Horsebean", 'Linseed', 'Meatmeal', 'Soybean', 'Sunflower'))

ggsave('chicksfeedtype.png', height = 6, width = 8, unit = 'in', dpi = 300)



#---------------------------------

data('sleep')

