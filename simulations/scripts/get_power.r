###########################################################################################
#########A function to calculate power from a dataframe with neutral, BS 1my and 3mya data#
#########Date of creation: 27.10/.2013	###################################################
#########Barbara Bitarello################
#########Last modified:28.10/2013#########
##########################################


get_power<-function(x, false_pos=0.05,nsims=1000){

	x[1,2]->a  #'neutral'
	x[1001,2]->b #'1 my'
	x[2001,2]->d #'3 my'
	table(is.na(x$NCV[1:1000]))[[1]]->n  #'how many actual simulations (neutral)
	false_pos*n->fp #number of false positives for a given false_pos rate
	1-false_pos->tn_rate #1-false_pos=specificity
	tn=tn_rate*n
	quantile(x$NCV[x$Pop==a], na.rm=T, probs=seq(0,1,0.01))[[(false_pos*100)+1]]->thr #find threshold value
	
	#power for 1 my
	table(is.na(x$NCV[1001:2000]))[[1]]->n1 #'how many actual simulations for 1 my
	table(x$NCV[x$Pop==b]<=thr)[[2]]->tp1  #how many below the thr
	tp1/n1->tp1r
	1-tp1r->fneg1
	#3mya
	table(is.na(x$NCV[2001:3000]))[[1]]->n3 #'how many actual simulations for 3 my
	table(x$NCV[x$Pop==d]<=thr)[[2]]->tp3 #how many below the thr
	tp3/n3->tp3r
	1-tp3r->fneg3

	l<-list(threshold_value=thr,false_pos=false_pos,n=n,n1=n1,tp1=tp1,tp1r=tp1r,fneg1=fneg1,n3=n3,tp3=tp3,tp3r=tp3r,fneg3=fneg3)
	return(l)
	}

##
