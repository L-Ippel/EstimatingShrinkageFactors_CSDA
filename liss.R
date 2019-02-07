###########################
##Liss panel data study
##predicting attrition 
##raw data file: Peter Lugtig 18-2-2014
##clean data
###########################

##read spss file:
library("foreign")
liss.dat		<- read.spss("Participation_panelmember_fieldworkperiod.sav", use.value.labels = F, to.data.frame = T)

#count how often respondents received a questionnaire: 
liss.dat$count	<- apply(liss.dat[,7:56], 1, FUN=function(x) sum(x!=-1,na.rm=T) )

#select the respondents that got at least one questionnaire
my.dat.wide		<- as.data.frame(liss.dat[which(liss.dat$count!=0),])

#in order to get some insight in how often respondents react
my.dat.wide$som	<- apply(my.dat.wide[,7:56], 1, FUN=function(x){sum(x==1)})
my.dat.wide$prob	<- my.dat.wide$som/my.dat.wide$count

#some respondents apprear to have dropped out before 11/2007, these respondents have te be removed from the data set: 
my.dat.wide$select<- 1
for(i in 1:nrow(my.dat.wide))
{
	if(my.dat.wide$stoplid[i]>200711|is.na(my.dat.wide$stoplid[i])) my.dat.wide$select[i] <- 1 
	else {my.dat.wide$select[i] <- 0}
}

my.dat.wide		<- my.dat.wide[my.dat.wide$select==1,]
my.dat.wide		<- my.dat.wide[my.dat.wide$startlid!=201111,]	#we only want the respondents that started before the last questionnaire
my.dat.wide$stoplid[which(my.dat.wide$stoplid==201107)] <- 201108 #201107 no one answered this questionnaire so this one will be deleted 

#for the analysis we want to take into account when someone is denoted as a drop out, 
#therefore we change the month labels by questionnaire numbers
nummering<-cbind(1:50, c(200711,200712,200801:200812,200901:200912,201001:201012,201101:201112))
my.dat.wide$stop	<- 999
for(i in 1:nrow(my.dat.wide))
{
	if(is.na(my.dat.wide$stoplid[i])){my.dat.wide$stop[i] <- NA}
	else {my.dat.wide$stop[i] <- nummering[which(my.dat.wide$stoplid[i]==nummering[,2]),1]}
}


my.dat.wide	<- my.dat.wide[!is.na(my.dat.wide[,2]),]	#some people dont have an id number (yet?)
my.dat.wide	<- my.dat.wide[,-51]				#july 2011 has 0 responses
my.dat.wide <- my.dat.wide[order(my.dat.wide[,2]),]	


write.table(my.dat.wide, "liss_wide.txt", sep=" ")
#the analysis is will be a replay of a data stream, for that purpose we change the data file from wide to long:
my.dat.long <- data.frame(id=numeric(0), response=numeric(0), 
			vragenlijst_nummer=numeric(0), prob = numeric(0), stop=numeric(0))
#7th variable is november 2007, start data set
my.dat.long <- rbind(my.dat.long,cbind(id=my.dat.wide[my.dat.wide[,7]!=-1,2], 
			response=my.dat.wide[my.dat.wide[,7]!=-1,7],
			vragenlijst_nummer=rep(1, length(my.dat.wide[my.dat.wide[,7]!=-1,7])),
			prob=my.dat.wide[my.dat.wide[,7]!=-1,58], stop=my.dat.wide[my.dat.wide[,7]!=-1,60]	))

#we randomly order the respondents within a questionnaire, order of the questionnaires is maintained:
set.seed(864531)
for(i in 1:48)	#50 questionnaires, the first is used to initiate the dataset, 1 questionnaire is dropped because there werent any responses
{
	nieuwe_lijst<-my.dat.wide[my.dat.wide[,7+i]!=-1,c(2,7+i,58,60)]	
	new.order<- sample(nrow(nieuwe_lijst)  )
	nieuwe_mensen<-nieuwe_lijst[new.order,]
	
	my.dat.long<-rbind(my.dat.long, cbind(id=nieuwe_mensen[,1], 
			response=nieuwe_mensen[,2],
			vragenlijst_nummer=rep(i+1, nrow(nieuwe_mensen)),
			prob=nieuwe_mensen[,3], stop=nieuwe_mensen[,4]))
}
my.dat.long	<- my.dat.long[-which(my.dat.long$stop<my.dat.long$vragenlijst_nummer),]

write.table(my.dat.long, "liss_long.txt", sep=" ")