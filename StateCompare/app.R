# A little shiny app to compare the log growth estimates of two segments of coronavirus case number data by US state
# Written by Stephen Proulx sproul@ucsb.edu
# data from https://covidtracking.com/about-tracker/

library(shiny)

library(tidyverse)
library(lubridate)
library(rstan)
options(mc.cores = parallel::detectCores())
library(bayesplot)


confirmed_sheet<-read.csv("http://covidtracking.com/api/states/daily.csv")  


confirmed_long <-  confirmed_sheet%>%
    mutate(date=date(dateChecked))%>% 
    mutate(delta.days=period_to_seconds(days(ymd(date) - ymd(20200122)))/(60*60*24)) %>% #calculate days since data began being collected
    rename(cases=positive)
    as_tibble()  


maxdays=max(confirmed_long$delta.days)



#make list of countries by most cases
state_list <- confirmed_long %>% filter(delta.days==maxdays) %>% arrange(desc(cases)) %>% select(state)
#countries in alphabetical order but only with more than 100 cases
state_list_alpha <- confirmed_long %>% filter(delta.days==maxdays,cases>25) %>% arrange(state) %>% select(state)

#mod <- readRDS("model.rds") 



# Define UI
ui <- fluidPage(
    
    # Application title
    titlePanel("Compare two states"),
    
    
    sidebarLayout(
        sidebarPanel(
            p("Day one is 01/22/2020."),
            selectInput("state1", "First state to compare:",
                        state_list,
                        selected = "CA"),
            selectInput("state2", "Seconed state to compare:",
                        state_list_alpha,
                        selected = "NY"),
            
            sliderInput("start1",
                        "Earliest day for state 1:",
                        min = 1,
                        max = maxdays,
                        value = 38),
            sliderInput("end1",
                        "Latest day for state 1:",
                        min = 1,
                        max = maxdays,
                        value = maxdays),
            sliderInput("start2",
                        "Earliest day for state 2:",
                        min = 1,
                        max = maxdays,
                        value = 45),
            sliderInput("end2",
                        "Latest day for state 2:",
                        min = 1,
                        max = maxdays,
                        value = maxdays),
            sliderInput("pstrength",
                        "Advanced: Strength of the prior that the two periods are the same",
                        round = -2, step = 0.01,
                        min = 0,
                        max = 2,
                        value = 0.5),
            actionButton("do", "Run Models")
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
            plotOutput("stateplot"),plotOutput("compare")
        )
    )
)



# Define server logic 
server <- function(input, output) {
    
    observeEvent(input$do, {
        tocomp = tibble(state = c(input$state1,input$state2) , startday = c(input$start1,input$start2) , 
                        endday = c(input$end1,input$end2)) %>%
            mutate(n=endday-startday)
        
        
        mydata.sub1 <- filter(confirmed_long, delta.days>=tocomp$startday[1],delta.days<=tocomp$endday[1],state == tocomp$state[1] , cases>25 ) %>% mutate(dataset=1)   
        mydata.sub2 <- filter(confirmed_long, delta.days>=tocomp$startday[2],delta.days<=tocomp$endday[2],state == tocomp$state[2] , cases>25 ) %>% mutate(dataset=2)
        
        mindays=c(min(mydata.sub1$delta.days),min(mydata.sub2$delta.days))
        
        mydata.sub1<-mydata.sub1 %>% mutate(days=delta.days-mindays[1]+1)
        mydata.sub2<-mydata.sub2 %>% mutate(days=delta.days-mindays[2]+1)
        
        maxdays=c(max(mydata.sub1$days),max(mydata.sub2$days))
        
        mydata.sub.total <-bind_rows(mydata.sub1,mydata.sub2) %>% group_by(dataset) %>% arrange(days, .by_group = TRUE) %>%
            rename(total_cases=cases) 
        
        total_datapoints=nrow(mydata.sub.total)
        
        stan_data <- c(mydata.sub.total[c("days","total_cases","dataset")],tocomp[c("state")],list(maxdays=maxdays,total_datapoints=total_datapoints , pstrength = input$pstrength)) 
        
        
        #    fit_poissCompare <- sampling(mod, 
        #                    data =stan_data, chains = 8,iter = 5000, seed = 2131231 )
        
        fit_poissCompare <- stan(file = 'model_compare2.stan', 
                                 data =stan_data, chains = 8,iter = 5000, seed = 2131231 )
        
        
        tmp<-as.data.frame(rstan::extract(fit_poissCompare,pars=c("mean_lambda","lambda_diff","lambda[1]","lambda[2]"))) %>% 
            rename("lambda1"="lambda.1.","lambda2"="lambda.2.")
        
        
        tmp1<-tmp %>%select(lambda1) %>%rename(mode=lambda1) %>%mutate(parameter="lambda1") %>% mutate(dtime=log(2)/(log(1+mode)))
        tmp2<-tmp %>%select(lambda2) %>%rename(mode=lambda2) %>%mutate(parameter="lambda2") %>% mutate(dtime=log(2)/(log(1+mode)))
        tmp<-rbind(tmp1,tmp2)
        
        
        lq=0.025
        uq=0.975
        d1<-density(tmp1$dtime)
        dd1 <- with(d1,data.frame(x,y)) %>% filter(y>0.01)
        qs1=quantile(ecdf(tmp1$dtime),prob=c(lq,0.5,uq))
        d2<-density(tmp2$dtime)
        dd2 <- with(d2,data.frame(x,y)) %>% filter(y>0.01)
        qs2=quantile(ecdf(tmp2$dtime),prob=c(lq,0.5,uq))
        
        maxdd1=max(dd1$y)
        dd1<-mutate(dd1,rely=y/maxdd1)
        maxdd2=max(dd2$y)
        dd2<-mutate(dd2,rely=y/maxdd2)
        
        upperlim=max(qs1[3],qs2[3])+0.5
        lowerlim=min(qs1[1],qs2[1])-0.5
        ddata<-bind_rows(mutate(tmp1,state=str_c(input$state1, " 1")),mutate(tmp2,state=str_c(input$state2, " 2"))) 
        
        output$compare <-  renderPlot({
        ggplot(data=ddata ) +
            geom_violin( aes(y=dtime,x=state ,  fill=state))+
            scale_y_continuous(limits=c(lowerlim,upperlim),name="Days until cases double")+
            scale_color_manual(breaks = c(str_c(input$state1, " 1"),str_c(input$state1, " 2")),values=c( "magenta","blue"))+
            annotate("text", x=c(0.8,1.8),
                     y=c(qs1[2]+0.25,qs2[2]+0.25), label= c(str_c(input$state1, " 1"),str_c(input$state2, " 2")),
                     size=3)+
            theme_bw()+
            theme( axis.text=element_text(family="Helvetica", size=8),text=element_text(family="Helvetica", size=12),legend.position = "none")
        
        })
        
        
        # output$compare <-  renderPlot({  
        #     boxdata=tibble(state=c(str_c(input$state1, " 1"), str_c(input$state1, " 2")), medians=c(qs1[2],qs2[2]) , lower=c(qs1[1],qs2[1]), upper=c(qs1[3],qs2[3]))
        #     
        #     ggplot(data=boxdata , aes(state,medians))+
        #         geom_point(sie=5) +
        #         geom_errorbar(aes(ymin=lower,ymax=upper),  width=0.5 )+
        #         scale_y_continuous(limits=c(0,10),name="Days until cases double")+
        #         annotate("text", x=c(1,2),
        #                  y=c(qs1[1]-.5,qs2[1]-.5), label= c(str_c(input$state1, " 1"),str_c(input$state2, " 2")),
        #                  size=5)+
        #          coord_fixed(ratio = 1/5) +
        #         theme_bw()+
        #         theme( axis.text=element_text(family="Helvetica", size=8),text=element_text(family="Helvetica", size=12))
        # })
        
        
        
        # output$compare <-  renderPlot({   
        #     ggplot(data=dd1 , aes(x,rely))+
        #         geom_line(data=dd1)+
        #         geom_ribbon(data=filter(dd1,x>qs1[[1]] & x<qs1[[3]]),aes(ymax=rely),ymin=0,
        #                     fill="#1B9E77",colour=NA,alpha=0.5)+
        #         geom_line(data=dd2)+
        #         geom_ribbon(data=filter(dd2,x>qs2[[1]] & x<qs2[[3]]),aes(ymax=rely),ymin=0,
        #                     fill="#D95F02",colour=NA,alpha=0.5)+
        #         scale_y_continuous(limits=c(0,1.5),name="Posterior Density")+
        #         scale_x_continuous(name="Days until cases double", limits=c(1, 8)) +
        #         annotate("text", x=c(qs1[[2]],qs2[[2]]),
        #                  y=1.2, label= c(str_c(input$state1, " 1"),str_c(input$state2, " 2")),
        #                  size=3)+
        #         theme_bw()+
        #         theme( axis.text=element_text(family="Helvetica", size=8),text=element_text(family="Helvetica", size=12))
        # })
        # 
        
        #    posterior <- as.matrix(fit_poissCompare)
        
        #    output$compare <-  renderPlot({mcmc_areas(posterior, pars = c("lambda_diff"),prob = 0.99) })
        
    })
    
    
    
    
    output$stateplot <- renderPlot({
        data1<-filter(confirmed_long,  state==input$state1,delta.days>(input$start1-1)  ,delta.days<(input$end1+1) ,cases>20) %>% 
            mutate(country=input$state1 ,dataset=1) %>% unite("state_set",country:dataset)
        data2<-filter(confirmed_long,  state==input$state2  ,delta.days>(input$start2-1)  ,delta.days<(input$end2+1) ,cases>20) %>% 
            mutate(country=input$state2 ,dataset=2) %>% unite("state_set",country:dataset)
        datatot<-bind_rows(data1,data2)
        
        ggplot(data=datatot , aes(x=delta.days,y=log(cases,base=10) , color=state_set,group=state_set)) +
            geom_point()+
            geom_smooth(method = "lm", formula = y ~ x) + 
            scale_y_continuous( limits=c(1,6),breaks=c(1,2,3,4,5), labels=c(10,100,1000,10000,100000))+
            scale_x_continuous( limits=c(20,maxdays))+
            labs( x="days since Jan 22" , y="cases", title="Log case numbers")
    })
    
   
    
    
    
}

# Run the application 
shinyApp(ui = ui, server = server)