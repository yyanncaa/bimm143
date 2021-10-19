#' ---
#' title: "Class 05"
#' author: "Byanca Lima Valenzuela 16688058"
#' date: "October 12, 2021"
#' ---
  
#Class 5: Data Visualization

# Lets start with a scatterplot
# before we can use it we need to load it up
library(ggplot2)

# every ggplot has a data + aes + geoms 
ggplot(data=cars) +
aes(x=speed,y=dist)+
geom_point()+
geom_smooth()

#change to a linear model
p<-ggplot(data=cars) +
aes(x=speed, y=dist)+
geom_point()+
geom_smooth(method="lm")
p

ggplot(cars) +
  aes(x=speed, y=dist)+
  geom_point() +
  labs(title="My  nice plot",
      x="Speed(MPH)",
      y="Stopping distance (ft)",
      subtitle="Your informativze subtitle text here",
      caption="Dataset: 'cars'")+
      geom_smooth(method="lm",se=FALSE)+
       theme_bw()





# Lets try a more complicated dataset of gene expression
# First read the dataset

url <- "https://bioboot.github.io/bimm143_S20/class-material/up_down_expression.txt"
genes <- read.delim(url)
head(genes)

#this code tells me how many genes
row(genes)
colnames(genes)
##[1]"Gene"   "Condition1" "Condition2" "State"
ncol(genes)
##[1]4
table(genes$State)
round(table(genes$State)/nrow(genes)*100,2)

p<-ggplot(genes)+
  aes(x=Condition1, y=Condition2, col=State)+
  geom_point()
p

p+scale_colour_manual(values=c("blue", "gray", "red"))+
  labs(title = "Gene Expression Chnages Upon Drug Treatment",
       x="Control(no drug)",
       y="Drug Treatment")